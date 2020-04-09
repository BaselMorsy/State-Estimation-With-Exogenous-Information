include("Ybus.jl")
import JuMP, Ipopt
using Parameters
using Distributions

module System_Definition
    using Parameters
    using DataFrames
    export System_Struct

    @with_kw mutable struct System_Struct
        y_bus
        b_line
        LineData
        BusData
        N_bus
        Sbase

        Operating_Cost = Float64[]
        state_error=Float64[]
        BusData_output = DataFrame(Bus = Int64[], V = Float64[], δ = Float64[],
         Pg = Float64[],Qg = Float64[],Pd = Float64[], Qd = Float64[])

        LineLoading = DataFrame(FromBus = Int64[],ToBus = Int64[]
            ,PL_1 = Float64[],PL_2 = Float64[],
            PLoss = Float64[],QL_1 = Float64[]
            ,QL_2 = Float64[],QLoss = Float64[])

        Line_Constraints= DataFrame(fbus = Union{Float64,Missing}[],
            tbus = Union{Float64,Missing}[], SL_Rating = Union{Float64,Missing}[])

        Measurements_PS = DataFrame(fbus=Int64[],tbus=Int64[],type = String[],Value=Float64[],Sigma = Float64[] )
        Measurements_Wind_true = []
        Measurements_Wind = []
        J=[]

        Estimated_Values = DataFrame(Bus=Int64[],V=Float64[],δ=Float64[])
        Estimated_wind_Speed = []
        Variable_Load_index = []
        Renewable_index = []

        Gen_Data =  DataFrame(bus = Int64[],
                Qmax = Union{Float64,Missing}[],
                Qmin = Union{Float64,Missing}[],
                Pmax = Union{Float64,Missing}[],
                Pmin = Union{Float64,Missing}[],
                C2 = Union{Float64,Missing}[],
                C1 = Union{Float64,Missing}[],
                C0 = Union{Float64,Missing}[])
        Pd_total = Float64[];

    end
end  # module

using Main.System_Definition;

function SystemAssembly(LineData,Sbase,BusData_sheet=nothing,
     GenData=nothing, Measurements_ps = nothing,
      Measurements_w = nothing)


    Lines = [];
    busData = [];
    YBUS  = [];
    bLine = [];
    genData = [];
    gen_constraints=[];
    S_rating_lines = [];
    N_bus =[];
    measurements_ps = [];
    measurements_w = [];


    if LineData != []
        Lines = DataFrame(CSV.read(LineData))
        N_bus = length(unique!(vcat(unique!(Vector(Lines[:fbus])),unique!(Vector(Lines[:tbus])))))
        YBUS,bLine = Y_Bus(Lines,N_bus)
        S_rating_lines = DataFrame(fbus = Lines[!,:fbus],tbus = Lines[!,:tbus]
        ,SL_Rating = Lines[!,:rate])
    end


    if BusData_sheet != nothing
        busData = DataFrame(CSV.read(BusData_sheet))
    end

    if GenData != nothing
        gendata = DataFrame(CSV.read(GenData))
        gen_data = DataFrame(bus = gendata[:bus],C2 = gendata[:c2],
            C1 = gendata[:c1],C0 = gendata[:c0],
            Pmax = gendata[:Pmax],
            Pmin = gendata[:Pmin],
            Qmax = gendata[:Qmax],
            Qmin = gendata[:Qmin]);
    end

    if Measurements_w != nothing
        measurements_w = DataFrame(CSV.read(Measurements_w));
        measurements_w[!,]
    end

    Sys = System_Struct(y_bus = YBUS,b_line = bLine, LineData = Lines,BusData =busData
        ,N_bus=N_bus,Sbase = Sbase)
    append!(Sys.Line_Constraints,S_rating_lines);
    Sys.Gen_Data=gen_data;
    Sys.Pd_total = sum(busData[:Pd]);
    Sys.Measurements_Wind_true = measurements_w;
    return Sys
end

function Batch_Solver(System ::System_Struct, Load_Index_Array, Load_TimeSeries,Renewable_Index,Renewable_TimeSeries)

    L_TS = DataFrame(CSV.read(Load_TimeSeries))./1000;
    R_TS = DataFrame(CSV.read(Renewable_TimeSeries))
    R_TS[!,2] = R_TS[!,2]./1000;
    N_Sim = size(L_TS,1)

    System.Variable_Load_index = Load_Index_Array;
    System.Renewable_index = Renewable_Index;

    System_array = Array{System_Struct}(undef,N_Sim)
    A = Array{System_Struct}(undef,N_Sim)

        [System_array[i] = System for i in 1:N_Sim]

        for i = 1:N_Sim

            k = 1;
            for j in Load_Index_Array

                System_array[i].BusData[j,3] = L_TS[i,k]
                System_array[i].BusData[j,4] = L_TS[i,k+1]
                k = k+2
            end
            m = 2;
            for l in Renewable_Index

                System_array[i].BusData[l,3] = -R_TS[i,m]
                System_array[i].Measurements_Wind_true = R_TS[!,1]
                m = m+2
            end

            A[i] = deepcopy(Solve_OPF(System_array[i]))
        end


        return A
end

function Solve_OPF(System ::System_Struct, method = "ACOPF", Contingency_Order = 1, Delta_P = 0)
    #Method can be "P-SCOPF", "C-SCOPF", or "ACOPF"
    #Contingency_Order is currently set to 1 -> N-1 Contingency
    G = real(System.y_bus);
    B = imag(System.y_bus);
    b = System.b_line;

    N = System.N_bus
    N_Lines = size(System.Line_Constraints,1)
    Nodes_set = 1:N;
    Gen_set = System.Gen_Data[!,:bus];

    Sbase = System.Sbase;

    BusData = System.BusData;
    GenData = System.Gen_Data;
    BranchData = System.Line_Constraints;
    S = zeros(N,N)

        for l in 1:N_Lines
            S[Int64(BranchData[l,:fbus]),Int64(BranchData[l,:tbus])] = BranchData[l,:SL_Rating]
            S[Int64(BranchData[l,:tbus]),Int64(BranchData[l,:fbus])] = BranchData[l,:SL_Rating]
        end

    if method == "ACOPF"


            m = JuMP.Model(JuMP.with_optimizer(Ipopt.Optimizer, print_level =0))

            # 2.1.Variables
            JuMP.@variable(m, BusData[i, :Vmin] ≤ v[i in Nodes_set] ≤ BusData[i, :Vmax])
            JuMP.@variable(m, -2*π ≤ δ[i in Nodes_set] ≤ 2*π)

            JuMP.@variable(m, GenData[(GenData[:bus].==g), :Pmin][1,1] ≤ p[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Pmax][1,1])
            JuMP.@variable(m, GenData[(GenData[:bus].==g), :Qmin][1,1] ≤ q[g in Gen_set] ≤ GenData[(GenData[:bus].==g), :Qmax][1,1])

            JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set])
            JuMP.@variable(m, qij[i = Nodes_set, j = Nodes_set])

            # 2.2. Constraints
            JuMP.@constraint(m, ReferenceAngle,
                (δ[1] ==  0.0))

            # ACTIVE POWER THROUGH LINE N-M
            JuMP.@NLconstraint(m, p_line[i in Nodes_set, j in Nodes_set],
                 (pij[i,j] ==  Sbase*(v[i]*v[j]*(G[i,j]*cos(δ[i]-δ[j])+B[i,j]*sin(δ[i]-δ[j])) -(v[i]^2)*G[i,j] )))

            # REACTIVE POWER THROUGH LINE N-M
            JuMP.@NLconstraint(m, q_line[i in Nodes_set, j in Nodes_set],
                 (qij[i,j] ==  Sbase*(v[i]*v[j]*(G[i,j]*sin(δ[i]-δ[j]) - B[i,j]*cos(δ[i]-δ[j])) +(v[i]^2)*(B[i,j]-b[i,j]) )));

            # ACTIVE NODAL BALANCE
            JuMP.@constraint(m, Pnodal[i in Nodes_set],
                sum(pij[i,j] for j = Nodes_set) == sum(p[g] for g in Gen_set if GenData[(GenData[:bus].==g),:bus][1,1]==i) - BusData[i,:Pd])

            # REACTIVE NODAL BALANCE
            JuMP.@constraint(m, Qnodal[i in Nodes_set],
                sum(qij[i,j] for j = Nodes_set) == sum(q[g] for g in Gen_set if GenData[(GenData[:bus].==g),:bus][1,1]==i) - BusData[i,:Qd])

            # LINE CAPACITY

            JuMP.@NLconstraint(m, Smax[i in Nodes_set, j in Nodes_set],
                pij[i,j]^2 + qij[i,j]^2 ≤ S[i,j]^2)

            #OBJECTIVE

            JuMP.@objective(m,Min,sum(GenData[(GenData[:bus].==g),:C2][1,1]*p[g]^2+GenData[(GenData[:bus].==g),:C1][1,1]*p[g]+GenData[(GenData[:bus].==g),:C0][1,1] for g in Gen_set))

            JuMP.optimize!(m)

            println("-------- AC OPTIMAL POWER FLOW -------")
                println("      Total Cost = "*string(round(JuMP.objective_value(m), digits =2))*" USD")
                println("-----------------------------------")
                println("Generation data:")
                for g in Gen_set
                    println("Gen "*string(g)*": P = "*string(round.(JuMP.value.(p)[g],digits =2))*" MW, Q = "*string(round.(JuMP.value.(q)[g],digits = 2))*" MVar")
                end

                println()
                Pg = JuMP.value.(p)
                Qg = JuMP.value.(q)
                Pij = JuMP.value.(pij)
                Qij = JuMP.value.(qij)
                V = JuMP.value.(v)
                δ = JuMP.value.(δ)
                Pd = System.BusData[!,3]
                Qd = System.BusData[!,4]

                System.Operating_Cost = JuMP.objective_value(m)
                Pg = [i in System.Gen_Data[!,1] ? Pg.data[findall(x -> x==i,System.Gen_Data[!,1])][1,1] : 0 for i in 1:System.N_bus]
                Qg = [i in System.Gen_Data[!,1] ? Qg.data[findall(x -> x==i,System.Gen_Data[!,1])][1,1] : 0 for i in 1:System.N_bus]
                df = DataFrame(Bus = 1:System.N_bus, V = V.data, δ = δ.data,
                 Pg = Pg ,Qg = Qg, Pd = Pd, Qd = Qd)
                System.BusData_output = df

                Pij_1 = [ Pij[System.Line_Constraints[i,1],System.Line_Constraints[i,2]] for i in 1:size(System.Line_Constraints,1) ]
                Pij_2 = [ Pij[System.Line_Constraints[i,2],System.Line_Constraints[i,1]] for i in 1:size(System.Line_Constraints,1) ]
                P_Loss = Pij_1 + Pij_2

                Qij_1 = [ Qij[System.Line_Constraints[i,1],System.Line_Constraints[i,2]] for i in 1:size(System.Line_Constraints,1) ]
                Qij_2 = [ Qij[System.Line_Constraints[i,2],System.Line_Constraints[i,1]] for i in 1:size(System.Line_Constraints,1) ]
                Q_Loss = Qij_1 + Qij_2
                df2 = DataFrame(FromBus =System.Line_Constraints[!,1] ,ToBus =System.Line_Constraints[!,2]
                    ,PL_1 = Pij_1,PL_2 = Pij_2,
                    PLoss = P_Loss,QL_1 = Qij_1
                    ,QL_2 = Qij_2,QLoss = Q_Loss)
                System.LineLoading = df2
                return System
    end

end

function Measurement_generator(System_array ::Array{System_Struct}, σ_P, σ_Q, σ_PL , σ_QL, σ_v, σ_X)

    for i in 1:length(System_array)
        df_PS = DataFrame(fbus = Int64[],tbus = Int64[],
            type = String[],value = Float64[],Sigma = Float64[])

        df_X = DataFrame(bus = Int64[],value = Float64[],Sigma=Float64[])

        bus_data = System_array[i].BusData_output
        line_data = System_array[i].LineLoading
        X_data = System_array[i].Measurements_Wind_true[i,:]

        #Nodal Measurements
        for j in 1:System_array[i].N_bus
            df_append_P=DataFrame(fbus = bus_data[j,:Bus],tbus = 0,type = "P_i"
                ,value = bus_data[j,:Pg]-bus_data[j,:Pd]+rand(Distributions.Normal(0,σ_P),1)[1,1]
                ,Sigma = σ_P)

            df_append_Q=DataFrame(fbus = bus_data[j,:Bus],tbus = 0,type = "Q_i"
                ,value = bus_data[j,:Qg]-bus_data[j,:Qd]+rand(Distributions.Normal(0,σ_Q),1)[1,1]
                ,Sigma = σ_Q)

            df_append_V=DataFrame(fbus = bus_data[j,:Bus],tbus = 0,type = "V_i"
                ,value = bus_data[j,:V] +rand(Distributions.Normal(0,σ_v),1)[1,1]
                ,Sigma = σ_v)

            append!(df_PS,df_append_V)
            append!(df_PS,df_append_P)
            append!(df_PS,df_append_Q)
        end

        #Line Measurements
        for k in 1:size(System_array[i].Line_Constraints,1)

            df_append_PL1=DataFrame(fbus = line_data[k,:FromBus]
                ,tbus = line_data[k,:ToBus],type = "PL_i"
                ,value = line_data[k,:PL_1] +rand(Distributions.Normal(0,σ_PL),1)[1,1]
                ,Sigma = σ_PL)
            df_append_PL2=DataFrame(fbus = line_data[k,:ToBus]
                ,tbus = line_data[k,:FromBus],type = "PL_i"
                ,value = line_data[k,:PL_2] +rand(Distributions.Normal(0,σ_PL),1)[1,1]
                ,Sigma = σ_PL)

            df_append_QL1=DataFrame(fbus = line_data[k,:FromBus]
                ,tbus = line_data[k,:ToBus],type = "QL_i"
                ,value = line_data[k,:QL_1] +rand(Distributions.Normal(0,σ_QL),1)[1,1]
                ,Sigma = σ_PL)
            df_append_QL2=DataFrame(fbus = line_data[k,:ToBus]
                ,tbus = line_data[k,:FromBus],type = "QL_i"
                ,value = line_data[k,:QL_2] +rand(Distributions.Normal(0,σ_QL),1)[1,1]
                ,Sigma = σ_PL)
            append!(df_PS,df_append_PL1)
            append!(df_PS,df_append_PL2)
            append!(df_PS,df_append_QL1)
            append!(df_PS,df_append_QL2)

        end

        #Renewable Measurements
        m = 1;
        for l in 1:length(System_array[i].Renewable_index)

            df_X_append = DataFrame(bus = System_array[i].Renewable_index[l],value = X_data[l,m] +rand(Distributions.Normal(0,σ_X),1)[1,1],
                            Sigma = σ_X)
            m = m + 2;
            append!(df_X,df_X_append)

        end
        System_array[i].Measurements_PS = deepcopy(df_PS)
        System_array[i].Measurements_Wind = deepcopy(df_X)
    end
end

function Solve_State_Estimation(System_array_Virtual,Include_Exogeneous ::Bool = false)

    if Include_Exogeneous
        k = (10^-6)*0.5*1.225*0.38*pi*0.25*56^2;
        println("SE with exogenous began!")
        for i in 1:length(System_array_Virtual)

            Nodes = 1:System_array_Virtual[i].N_bus
            NBuses = System_array_Virtual[i].N_bus
            Buses = Nodes
            b = System_array_Virtual[i].b_line
            GL = real(System_array_Virtual[i].y_bus)
            BL = imag(System_array_Virtual[i].y_bus)

            Pi_M,Qi_M,PL_M,QL_M,Vi_M,Vw_M = extract_measurements_PS(System_array_Virtual[i])

            StateSt = JuMP.Model(JuMP.with_optimizer( Ipopt.Optimizer,print_level =1))

        # Variables
        JuMP.@variable(StateSt, P[1:NBuses])
        JuMP.@variable(StateSt, Q[1:NBuses])


        JuMP.@variable(StateSt, PL[1:NBuses, 1:NBuses])
        JuMP.@variable(StateSt, QL[1:NBuses, 1:NBuses])

        JuMP.@variable(StateSt, V[1:NBuses] >= 0.0, start=1.0)
        JuMP.@variable(StateSt, -1*pi ≤ δ[1:NBuses] ≤ pi, start=0.0)

        JuMP.@variable(StateSt,V_w >= 0)

        JuMP.@variable(StateSt, TotalCost>=0)

        # Constraints
        JuMP.@constraint(StateSt, ReferenceAngle,
            (δ[1] ==  0.0));

        JuMP.@NLconstraint(StateSt, Pinjected[i in Buses],
            (P[i] ==  V[i]*sum( V[j]*(GL[i,j]*cos(δ[i]-δ[j]) + BL[i,j]*sin(δ[i]-δ[j])) for j in Buses)));

        JuMP.@NLconstraint(StateSt, Qinjected[i in Buses],
            (Q[i] ==  V[i]*sum( V[j]*(GL[i,j]*sin(δ[i]-δ[j]) - BL[i,j]*cos(δ[i]-δ[j])) for j in Buses)));

        JuMP.@NLconstraint(StateSt, Pline[i in Buses, j in Buses],
            (PL[i,j] ==  V[i]*V[j]*(GL[i,j]*cos(δ[i]-δ[j])+BL[i,j]*sin(δ[i]-δ[j])) -V[i]^2*GL[i,j] ));

        JuMP.@NLconstraint(StateSt, Qline[i in Buses, j in Buses],
            (QL[i,j] ==  V[i]*V[j]*(GL[i,j]*sin(δ[i]-δ[j]) - BL[i,j]*cos(δ[i]-δ[j])) +V[i]^2*(BL[i,j] - b[i,j]) ));

        # Constraints to include exogenous parameters
        JuMP.@NLconstraint(StateSt,Vwind,V_w^3 - (sum(PL[5,j] for j in [2,4])/k) <= ((sqrt(1/Vw_M[5,:W]))^(1/3)+sqrt(1/Pi_M[5, :W])))
        JuMP.@NLconstraint(StateSt,Vwind2,(sum(PL[5,j] for j in [2,4])/k) - V_w^3 <= ((sqrt(1/Vw_M[5,:W]))^(1/3)+sqrt(1/Pi_M[5, :W])))
        #JuMP.@NLconstraint(StateSt,Vwind,V_w^3 == (sum(PL[5,j] for j in [2,4])/k))

        JuMP.@NLconstraint(StateSt, ObjectiveFunction,
             (TotalCost ==
             sum(Vi_M[i, :W]*(V[i] - Vi_M[i, :Value])^2 for i in Buses if ((i,:Value) in keys(Vi_M)) )
             + sum(Pi_M[i, :W]*(P[i] - Pi_M[i, :Value])^2 for i in Buses if ((i,:Value) in keys(Pi_M)) )
             + sum(Qi_M[i, :W]*(Q[i] - Qi_M[i, :Value])^2 for i in Buses if ((i,:Value) in keys(Qi_M)) )
             + sum(QL_M[i, j, :W]*(QL[i, j] - QL_M[i, j, :Value])^2 for i in Buses, j in Buses if ((i,j,:Value) in keys(QL_M))  )
             + sum(PL_M[i, j, :W]*(PL[i, j] - PL_M[i, j, :Value])^2 for i in Buses, j in Buses if ((i,j,:Value) in keys(PL_M)) )
             + Vw_M[5,:W]*(V_w-Vw_M[5,:value])^2));



    #  Objective
    JuMP.@NLobjective(StateSt, Min, TotalCost);

    JuMP.optimize!(StateSt)

    Cost = JuMP.objective_value(StateSt)

     V = JuMP.value.(V)
     δ = JuMP.value.(δ)
     P = JuMP.value.(P)
     Q = JuMP.value.(Q)
    pij = JuMP.value.(PL)
    qij = JuMP.value.(QL)
    Vw = JuMP.value(V_w)

    println(Cost)
    System_array_Virtual[i].J = Cost;
    System_array_Virtual[i].Estimated_Values = DataFrame(Bus=1:System_array_Virtual[i].N_bus,V=V,δ=δ);
    System_array_Virtual[i].Estimated_wind_Speed = Vw
        end
        println("SE with exogenous ended!")
    else
        println("SE without exogenous began!")
        for i in 1:length(System_array_Virtual)
            Nodes = 1:System_array_Virtual[i].N_bus
            NBuses = System_array_Virtual[i].N_bus
            Buses = Nodes
            b = System_array_Virtual[i].b_line
            GL = real(System_array_Virtual[i].y_bus)
            BL = imag(System_array_Virtual[i].y_bus)

            Pi_M,Qi_M,PL_M,QL_M,Vi_M,Vw_M = extract_measurements_PS(System_array_Virtual[i])

            StateSt = JuMP.Model(JuMP.with_optimizer( Ipopt.Optimizer, print_level =1))

        # Variables
        JuMP.@variable(StateSt, P[1:NBuses])
        JuMP.@variable(StateSt, Q[1:NBuses])


        JuMP.@variable(StateSt, PL[1:NBuses, 1:NBuses])
        JuMP.@variable(StateSt, QL[1:NBuses, 1:NBuses])

        JuMP.@variable(StateSt, V[1:NBuses] >= 0.0, start=1.0)
        JuMP.@variable(StateSt, -1*pi ≤ δ[1:NBuses] ≤ pi, start=0.0)

        JuMP.@variable(StateSt, TotalCost>=0)

        # Constraints
        JuMP.@constraint(StateSt, ReferenceAngle,
            (δ[1] ==  0.0));

        JuMP.@NLconstraint(StateSt, Pinjected[i in Buses],
            (P[i] ==  V[i]*sum( V[j]*(GL[i,j]*cos(δ[i]-δ[j]) + BL[i,j]*sin(δ[i]-δ[j])) for j in Buses)));

        JuMP.@NLconstraint(StateSt, Qinjected[i in Buses],
            (Q[i] ==  V[i]*sum( V[j]*(GL[i,j]*sin(δ[i]-δ[j]) - BL[i,j]*cos(δ[i]-δ[j])) for j in Buses)));

        JuMP.@NLconstraint(StateSt, Pline[i in Buses, j in Buses],
            (PL[i,j] ==  V[i]*V[j]*(GL[i,j]*cos(δ[i]-δ[j])+BL[i,j]*sin(δ[i]-δ[j])) -V[i]^2*GL[i,j] ));

        JuMP.@NLconstraint(StateSt, Qline[i in Buses, j in Buses],
            (QL[i,j] ==  V[i]*V[j]*(GL[i,j]*sin(δ[i]-δ[j]) - BL[i,j]*cos(δ[i]-δ[j])) +V[i]^2*(BL[i,j] - b[i,j]) ));

        JuMP.@NLconstraint(StateSt, ObjectiveFunction,
             (TotalCost ==
             sum(Vi_M[i, :W]*(V[i] - Vi_M[i, :Value])^2 for i in Buses if ((i,:Value) in keys(Vi_M)) )
             + sum(Pi_M[i, :W]*(P[i] - Pi_M[i, :Value])^2 for i in Buses if ((i,:Value) in keys(Pi_M)) )
             + sum(Qi_M[i, :W]*(Q[i] - Qi_M[i, :Value])^2 for i in Buses if ((i,:Value) in keys(Qi_M)) )
             + sum(QL_M[i, j, :W]*(QL[i, j] - QL_M[i, j, :Value])^2 for i in Buses, j in Buses if ((i,j,:Value) in keys(QL_M))  )
             + sum(PL_M[i, j, :W]*(PL[i, j] - PL_M[i, j, :Value])^2 for i in Buses, j in Buses if ((i,j,:Value) in keys(PL_M)) )
             ));

    # Objective
    JuMP.@NLobjective(StateSt, Min, TotalCost);

    JuMP.optimize!(StateSt)

    Cost = JuMP.objective_value(StateSt)

     V = JuMP.value.(V)
     δ = JuMP.value.(δ)
     P = JuMP.value.(P)
     Q = JuMP.value.(Q)
    pij = JuMP.value.(PL)
    qij = JuMP.value.(QL)

    println(Cost)
    System_array_Virtual[i].J = Cost;
    System_array_Virtual[i].Estimated_Values = DataFrame(Bus=1:System_array_Virtual[i].N_bus,V=V,δ=δ);
        end
        println("SE without exogenous ended!")
    end
end

function extract_measurements_PS(System ::System_Struct)

    Measurements = System.Measurements_PS
    Measurements_W = System.Measurements_Wind

    pij_M = Dict{}()

    qij_M = Dict{}()

    Pi_M = Dict{}()

    Qi_M = Dict{}()

    Vi_M = Dict{}()

    Vw_M = Dict{}()

    for i in 1:size(Measurements,1)
        Type = Measurements[i,:type]
        FromBus = Measurements[i,:fbus]
        ToBus = Measurements[i,:tbus]
        Value = Measurements[i,:value]
        W = Measurements[i,:Sigma]

        if Type == "PL_i"
            pij_M[FromBus, ToBus, :Value] = Value
            pij_M[FromBus, ToBus, :W] = 1/W^2

        elseif Type == "QL_i"
            qij_M[FromBus, ToBus, :Value] = Value
            qij_M[FromBus, ToBus, :W] = 1/W^2

        elseif Type == "P_i"
            Pi_M[FromBus, :Value] = Value
            Pi_M[FromBus, :W] = 1/W^2

        elseif Type == "Q_i"
            Qi_M[FromBus, :Value] = Value
            Qi_M[FromBus, :W] = 1/W^2

        elseif Type == "V_i"
            Vi_M[FromBus, :Value] = Value
            Vi_M[FromBus, :W] = 1/W^2

        end
    end

    for i in 1:size(Measurements_W,1)
        Value = Measurements_W[i,:value]
        W = Measurements_W[i,:Sigma]
        Bus = Measurements_W[i,:bus]

        Vw_M[Bus,:value] = Value
        Vw_M[Bus,:W] = 1/W^2

    end

    return Pi_M,Qi_M,pij_M,qij_M,Vi_M,Vw_M
end

function calculate_residuals(System ::System_Struct)

    N_bus = System.N_bus
    V = System.Estimated_Values.V
    δ = System.Estimated_Values.δ
    G = real(System.y_bus)
    B = imag(System.y_bus)
    LineData = System.LineData

    LineIndices = [(LineData[i,1:2][1],LineData[i,1:2][2]) for i in 1:size(LineData,1)]

    Pi_M,Qi_M,pij_M,qij_M,Vi_M,Vw_M=extract_measurements_PS(System)
    Vi_M = [Vi_M[(i, :Value)] for i in 1:N_bus]
    Pi_M = [Pi_M[(i, :Value)] for i in 1:N_bus]
    Qi_M = [Qi_M[(i, :Value)] for i in 1:N_bus]
    pij_M_1 = [pij_M[(i[1], i[2], :Value)] for i in LineIndices]
    pij_M_2 = [pij_M[(i[2], i[1], :Value)] for i in LineIndices]
    qij_M_1 = [qij_M[(i[1], i[2], :Value)] for i in LineIndices]
    qij_M_2 = [qij_M[(i[2], i[1], :Value)] for i in LineIndices]
    Z = vcat(Vi_M,Pi_M,Qi_M,pij_M_1,pij_M_2,qij_M_1,qij_M_2,Vw_M[(5, :value)])

    Pi = [get_P_i(V,δ,G,B,i) for i in 1:N_bus]
    Qi = [get_Q_i(V,δ,G,B,i) for i in 1:N_bus]
    PL_1 = [get_P_L(V,δ,G,B,i[1],i[2]) for i in LineIndices]
    PL_2 = [get_P_L(V,δ,G,B,i[2],i[1]) for i in LineIndices]
    QL_1 = [get_Q_L(V,δ,G,B,LineData,i[1],i[2]) for i in LineIndices]
    QL_2 = [get_Q_L(V,δ,G,B,LineData,i[2],i[1]) for i in LineIndices]
    h = vcat(V,Pi,Qi,PL_1,PL_2,QL_1,QL_2,System.Estimated_wind_Speed)

    r = Z-h

    return r, Z, h
end

function calculate_cov_matrix(System ::System_Struct, IncludeExogenous = false)

    V = System.Estimated_Values.V
    δ = System.Estimated_Values.δ
    G = real(System.y_bus)
    B = imag(System.y_bus)
    LineData = System.LineData
    N_bus = System.N_bus
    LineIndices = [(LineData[i,1:2][1],LineData[i,1:2][2]) for i in 1:size(LineData,1)]

    Measurements = System.Measurements_PS

    Vi_M = @linq Measurements |> where(:type .== "V_i") |>
        select(Bus=:fbus, :value, :Sigma)

    Pi_M = @linq Measurements |> where(:type .== "P_i") |>
        select(Bus=:fbus, :value, :Sigma)

    Qi_M = @linq Measurements |> where(:type .== "Q_i") |>
        select(Bus=:fbus, :value, :Sigma)

    PL_M = @linq Measurements |> where(:type .== "PL_i") |>
        select(:fbus, :tbus, :value, :Sigma)

    QL_M = @linq Measurements |> where(:type .== "QL_i") |>
        select(:fbus, :tbus, :value, :Sigma)

    M_Sizes = [size(Vi_M,1) size(Pi_M,1) size(Qi_M,1) size(PL_M,1) size(QL_M,1)]
    M = (Vi_M, Pi_M, Qi_M, PL_M, QL_M)

    J,h = Generate_Jacobian(System,V,δ,G,B,LineData,M_Sizes,M,IncludeExogenous)


    Pi_M,Qi_M,pij_M,qij_M,Vi_M,Vw_M=extract_measurements_PS(System)
    Vi_M_W = [Vi_M[(i, :W)] for i in 1:N_bus]
    Pi_M_W = [Pi_M[(i, :W)] for i in 1:N_bus]
    Qi_M_W = [Qi_M[(i, :W)] for i in 1:N_bus]
    pij_M_W_1 = [pij_M[(i[1], i[2], :W)] for i in LineIndices]
    pij_M_W_2 = [pij_M[(i[2], i[1], :W)] for i in LineIndices]
    qij_M_W_1 = [qij_M[(i[1], i[2], :W)] for i in LineIndices]
    qij_M_W_2 = [qij_M[(i[2], i[1], :W)] for i in LineIndices]
    M_vec = vcat(Vi_M_W,Pi_M_W,Qi_M_W,pij_M_W_1,pij_M_W_2,qij_M_W_1,qij_M_W_2)

    if IncludeExogenous
        M_vec = vcat(M_vec,Vw_M[(5, :W)])
        W = Diagonal(M_vec)
    else
        W = Diagonal(M_vec)
    end

    Ω = inv(W) - J*inv(transpose(J)*W*J)*transpose(J)

    return Ω,J
end

# Helper functions
function get_P_i(V,δ,G,B,i)

    Pi = V[i];
    S = 0;

    for j in 1:length(V)
        S = S + V[j]*(G[i,j]*cos(δ[i]-δ[j]) + B[i,j]*sin(δ[i]-δ[j]));
    end
    Pi = Pi*S;
    return Pi
end

function get_P_L(V,δ,G,B,i,j)

    if i != j
        P_L = V[i]*V[j]*(G[i,j]*cos(δ[i]-δ[j])+B[i,j]*sin(δ[i]-δ[j])) - (V[i]^2)*G[i,j];
    else
        P_L = [];
    end
     return P_L
 end

function get_Q_i(V,δ,G,B,i)
     Qi = V[i];
     S = 0;

     for j in 1:length(V)
         S = S + V[j]*(G[i,j]*sin(δ[i]-δ[j]) - B[i,j]*cos(δ[i]-δ[j]));
     end
     Qi = Qi*S;
     return Qi
end

function get_Q_L(V,δ,G,B,LineData,i,j)

        if i != j
                b1 = @linq LineData |> where(:fbus .== i, :tbus .== j) |>
                select(:b);

                if ~isempty(b1)
                        b1 = b1[1,1];
                else
                        b1 = 0;
                end

                b2 = @linq LineData |> where(:fbus .== j, :tbus .== i) |>
                select(:b);

                if ~isempty(b2)
                        b2 = b1[1,1];
                else
                        b2 = 0;
                end

                b = b1 + b2;

                Q_L = P_L = V[i]*V[j]*(G[i,j]*sin(δ[i]-δ[j])-B[i,j]*cos(δ[i]-δ[j])) + (V[i]^2)*(B[i,j] - b);
         else
                 Q_L = [];
         end

 return Q_L
end

function Generate_Jacobian(System,V,δ,G,B,LineData,M_Sizes,M,IncludeExogenous = false,wind_bus = 5)
    a = M_Sizes[1];
    b = M_Sizes[2];
    c = M_Sizes[3];
    d = M_Sizes[4];
    e = M_Sizes[5];
    N = length(V);

    #Voltage part of the jacobian
    #V_V
    V_V = zeros(a,N)*1.0;
    if !isempty(M[1])
        for i in 1:a
            b_i = M[1][:Bus][i];
            for j in 1:N
                if j == b_i
                    V_V[i,j] = 1;
                else
                    V_V[i,j] = 0;
                end
            end
        end
    end

    V_δ = zeros(a,N-1)*1.0; #V_δ

    #P_V part
    P_V = zeros(b,N)*1.0;
    if !isempty(M[2])
        for i in 1:b
            b_i = M[2][:Bus][i];
            for j in 1:N
                if b_i == j #Pi_Vi
                    P_V[i,j] = V[b_i]*G[b_i,j];
                    for k in 1:N
                        P_V[i,j] += V[k]*(G[b_i,k]*cos(δ[b_i]-δ[k])+B[b_i,k]*sin(δ[b_i]-δ[j]));
                    end
                else #Pi_Vj
                    P_V[i,j] = V[b_i]*(G[b_i,j]*cos(δ[b_i]-δ[j])+B[b_i,j]*sin(δ[b_i]-δ[j]));
                end
            end
        end
    end

    #P_δ part
    P_δ = zeros(b,N-1)*1.0;
    if !isempty(M[2])
        for i in 1:b
            b_i = M[2][:Bus][i];
            for j in 2:N
                if b_i == j #Pi_δi
                    P_δ[i,j-1] = -(V[b_i]^2)*B[b_i,j];
                    for k in 1:N
                        P_δ[i,j-1] += V[b_i]*V[k]*(-G[b_i,k]*sin(δ[b_i]-δ[k])+B[b_i,k]*cos(δ[b_i]-δ[j]));
                    end
                else #Pi_δj
                    P_δ[i,j-1] = V[j]*V[b_i]*(G[b_i,j]*sin(δ[b_i]-δ[j])-B[b_i,j]*cos(δ[b_i]-δ[j]));
                end
            end
        end
    end

    #Q_V
    Q_V = zeros(c,N)*1.0;
    if !isempty(M[3])
        for i in 1:c
            b_i = M[3][:Bus][i]
            for j in 1:N
                if b_i == j #Qi_Vi
                    Q_V[i,j] = -V[b_i]*B[b_i,j];
                    for k in 1:N
                        Q_V[i,j] += V[k]*(G[b_i,k]*sin(δ[b_i]-δ[k])-B[b_i,k]*cos(δ[b_i]-δ[j]));
                    end

                else #Qi_Vj
                    Q_V[i,j] = V[b_i]*(G[b_i,j]*sin(δ[b_i]-δ[j])-B[b_i,j]*cos(δ[b_i]-δ[j]));
                end
            end
        end
    end

    #Q_δ
    Q_δ = zeros(c,N-1)*1.0;
    if !isempty(M[3])
        for i in 1:c
            b_i = M[3][:Bus][i]
            for j in 2:N
                if i == j #Qi_δi
                    Q_δ[i,j-1] = -(V[b_i]^2)*G[b_i,j];
                    for k in 1:N
                        Q_δ[i,j-1] += V[b_i]*V[k]*(G[b_i,k]*cos(δ[b_i]-δ[k])+B[b_i,k]*sin(δ[b_i]-δ[j]));
                    end
                else #Qi_δj
                    Q_δ[i,j-1] = -V[j]*V[b_i]*(G[b_i,j]*cos(δ[b_i]-δ[j])+B[b_i,j]*sin(δ[b_i]-δ[j]));
                end
            end
        end
    end

    #PL_V
    PL_V = zeros(d,N)*1.0;
    if !isempty(M[4])
        for i in 1:d
            b1_i = M[4][:fbus][i];
            b2_j = M[4][:tbus][i];
            for j in 1:N
                if j == b1_i #Pij_Vi
                    PL_V[i,j] = V[b2_j]*(G[b1_i,b2_j]*cos(δ[b1_i]-δ[b2_j])
                    +B[b1_i,b2_j]*sin(δ[b1_i]-δ[b2_j]))-2*G[b1_i,b2_j]*V[b1_i];
                elseif j == b2_j #Pij_Vj
                    PL_V[i,j] = V[b1_i]*(G[b1_i,b2_j]*cos(δ[b1_i]-δ[b2_j])
                    +B[b1_i,b2_j]*sin(δ[b1_i]-δ[b2_j]))
                else
                    PL_V[i,j] = 0.0;
                end
            end
        end
    end

    #PL_δ
    PL_δ = zeros(d,N-1)*1.0;
    if !isempty(M[4])
        for i in 1:d
            b1_i = M[4][:fbus][i];
            b2_j = M[4][:tbus][i];
            for j in 2:N
                if j == b1_i #Pij_δi
                    PL_δ[i,j-1] = V[b1_i]*V[b2_j]*(-G[b1_i,b2_j]*sin(δ[b1_i]-δ[b2_j])
                    +B[b1_i,b2_j]*cos(δ[b1_i]-δ[b2_j]));
                elseif j == b2_j #Pij_δj
                    PL_δ[i,j-1] = V[b1_i]*V[b2_j]*(G[b1_i,b2_j]*sin(δ[b1_i]-δ[b2_j])
                    -B[b1_i,b2_j]*cos(δ[b1_i]-δ[b2_j]));
                else
                    PL_δ[i,j-1] = 0.0;
                end
            end
        end
    end

    #QL_V
    QL_V = zeros(e,N)*1.0
    if !isempty(M[5])
        for i in 1:c
            b1_i = M[5][:fbus][i];
            b2_j = M[5][:tbus][i];
            for j in 1:N
                if j == b1_i #Qij_Vi

                    b = @linq LineData |> where(:fbus .== b1_i,:tbus .== b2_j) |> select(:b);
                    if size(b,1) == 0
                        b = @linq LineData |> where(:fbus .== b2_j,:tbus .== b1_i) |> select(:b);
                    end

                    QL_V[i,j] = V[b2_j]*(G[b1_i,b2_j]*sin(δ[b1_i]-δ[b2_j])-B[b1_i,b2_j]*cos(δ[b1_i]-δ[b2_j]))+2*(B[b1_i,b2_j]-b[1,1])*V[b1_i];
                elseif j == b2_j #Qij_Vj
                    QL_V[i,j] = V[b1_i]*(G[b1_i,b2_j]*sin(δ[b1_i]-δ[b2_j])-B[b1_i,b2_j]*cos(δ[b1_i]-δ[b2_j]));
                else
                    QL_V[i,j] = 0.0;
                end
            end
        end
    end

    #QL_δ
    QL_δ = zeros(e,N-1)*1.0
    if !isempty(M[5])
        for i in 1:d
            b1_i = M[5][:fbus][i];
            b2_j = M[5][:tbus][i];
            for j in 2:N
                if j == b1_i #Pij_δi
                    QL_δ[i,j-1] = V[b1_i]*V[b2_j]*(G[b1_i,b2_j]*cos(δ[b1_i]-δ[b2_j])
                    +B[b1_i,b2_j]*sin(δ[b1_i]-δ[b2_j]));
                elseif j == b2_j #Pij_δj
                    QL_δ[i,j-1] = -V[b1_i]*V[b2_j]*(G[b1_i,b2_j]*cos(δ[b1_i]-δ[b2_j])
                    +B[b1_i,b2_j]*sin(δ[b1_i]-δ[b2_j]));;
                else
                    QL_δ[i,j-1] = 0.0;
                end
            end
        end
    end

    J=[V_V V_δ; P_V P_δ; Q_V Q_δ; PL_V PL_δ; QL_V QL_δ];

    if IncludeExogenous
        k = (10^-6)*0.5*1.225*0.38*pi*0.25*56^2;
        Vw_row = J[a+wind_bus,:]./(3*k*System.Estimated_wind_Speed^2)
        J = vcat(J,transpose(Vw_row))
    end

    #------------------------------------------------------#
    #---------------------Generate_h(X)--------------------#

    V_i = zeros(a,1)*1.0;

    #V_i
    for i in 1:a
        bus = M[1][:Bus][i];
        V_i[i] = V[bus];
    end

    P_i = reshape( !isempty(M[2]) ? [get_P_i(V,δ,G,B,i) for i in M[2][:Bus]] : [] ,:,1);
    Q_i = reshape( !isempty(M[3]) ? [get_Q_i(V,δ,G,B,i) for i in M[3][:Bus]] : [] ,: , 1);

    PL = reshape( !isempty(M[4]) ? [get_P_L(V,δ,G,B,i,j) for i in M[4][:fbus] , j in M[4][:tbus]] : [],:,1);
    QL = reshape( !isempty(M[5]) ? [get_Q_L(V,δ,G,B,LineData,i,j) for i in M[5][:fbus] , j in M[5][:tbus]] : [],:,1);

    h = vcat(V_i, P_i, Q_i, PL, QL)
    if IncludeExogenous
        h = vcat(h,System.Estimated_wind_Speed)
    end

    return J,h
end
#---------------------------------------#

function save_measurements(System,k)

    CSV.write(string("Measurements_",string(k),".csv"),System[k].Measurements_PS)
end
