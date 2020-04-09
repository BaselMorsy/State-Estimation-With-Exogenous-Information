using LinearAlgebra,DataFrames,DataFramesMeta,CSV;
import JuMP, Ipopt
import Distributions, Random
using Plots
using Statistics
using StatsPlots
using StatsFuns
using StatsBase
using Dates

include("SystemAssembly.jl")
include("report.jl")
Sbase = 1;

#Instantiate System
System = SystemAssembly("branch.csv",Sbase,"bus.csv","gen.csv")

#Solve OPF for all time instance
System_array = Batch_Solver(System,[2,3,4],"DemandPower.csv",[5],"wind.csv");


#Generate Measurements
α_P = 0.01
σ_Q = 0.01
σ_PL = 0.01
σ_QL =0.01
σ_V = 0.01
σ_X = 0.01

n=1
        E1_array=Vector{Float64}(undef,n)
        E2_array=Vector{Float64}(undef,n)

        E1_cheating_array=Vector{Float64}(undef,n)
        E2_cheating_array=Vector{Float64}(undef,n)

        MSE_V_array=Vector{Any}(undef,n)
        MSE_δ_array=Vector{Any}(undef,n)

        MSE_V_X_array=Vector{Any}(undef,n)
        MSE_δ_X_array=Vector{Any}(undef,n)

        MSE_V_cheating_array=Vector{Any}(undef,n)
        MSE_δ_cheating_array=Vector{Any}(undef,n)

        MSE_V_X_cheating_array=Vector{Any}(undef,n)
        MSE_δ_X_cheating_array=Vector{Any}(undef,n)

        J_arr = Vector{Any}(undef,n)
        J_X_arr = Vector{Any}(undef,n)

System_array_Virtual = deepcopy(System_array)
for i in 1:n

        #Generate Measurements
        Measurement_generator(System_array,α_P,σ_Q,σ_PL,σ_QL,σ_V,σ_X)

        #Solve State Estimation without Exogenous Parameters
        Solve_State_Estimation(System_array)
        J = [vcat(System_array[i].J) for i in 1: length(System_array)]

        #Solve State Estimation with Exogenous Parameters
        Solve_State_Estimation(System_array_Virtual,true)
        J_X = [vcat(System_array_Virtual[i].J) for i in 1: length(System_array_Virtual)]


        #Calculate squared errors for voltages
        SE_V_X = [vcat((System_array_Virtual[i].BusData_output[:V]-System_array_Virtual[i].Estimated_Values[:V]).^2) for i in 1:length(System_array_Virtual)]
        SE_V = [vcat((System_array_Virtual[i].BusData_output[:V]-System_array[i].Estimated_Values[:V]).^2) for i in 1:length(System_array_Virtual)]
        t = 1:length(System_array_Virtual)

        #Calulation of MSE for Voltages
        MSE_V = mean(SE_V[1:144])
        MSE_V_X = mean(SE_V_X[1:144])

        #Calculate Sum of MSE for Voltages
        SMSE_V = sum(MSE_V)
        SMSE_V_X = sum(MSE_V_X)

        #Calulation of MSE for Voltages in case of cheating
        MSE_V_cheating = mean(SE_V[144:192])
        MSE_V_X_cheating= mean(SE_V_X[144:192])

        #Calculate Sum of MSE for Voltages in case of cheating
        SMSE_V_cheating = sum(MSE_V_cheating)
        SMSE_V_X_cheating = sum(MSE_V_X_cheating)


        #Calculate squared errors for phase angles
        SE_δ_X = [vcat((System_array_Virtual[i].BusData_output[:δ]-System_array_Virtual[i].Estimated_Values[:δ]).^2) for i in 1:length(System_array_Virtual)]
        SE_δ = [vcat((System_array_Virtual[i].BusData_output[:δ]-System_array[i].Estimated_Values[:δ]).^2) for i in 1:length(System_array_Virtual)]


        #Calulation of MSE for phase angles
        MSE_δ = mean(SE_δ[1:144])
        MSE_δ_X = mean(SE_δ_X[1:144])

        #Calculate Sum of MSE for phase angles
        SMSE_δ = sum(MSE_δ)
        SMSE_δ_X = sum(MSE_δ_X)

        #Calulation of MSE for phase angles in case of cheating
        MSE_δ_cheating = mean(SE_δ[144:192])
        MSE_δ_X_cheating= mean(SE_δ_X[144:192])

        #Calculate Sum of MSE for phase angles in case of cheating
        SMSE_δ_cheating = sum(MSE_δ_cheating)
        SMSE_δ_X_cheating = sum(MSE_δ_X_cheating)

        #Calculation of error index for SE without and with exogenous
        E1=SMSE_V+SMSE_δ
        E2=SMSE_V_X+SMSE_δ_X

        #Calculation of error index for SE without and with exogenous in case of cheating
        E1_cheating = SMSE_V_cheating+SMSE_δ_cheating
        E2_cheating = SMSE_V_X_cheating+SMSE_δ_X_cheating


        #assignments:
        E1_array[i]=deepcopy(E1)
        E2_array[i]=deepcopy(E2)
        E1_cheating_array[i]=deepcopy(E1_cheating)
        E2_cheating_array[i]=deepcopy(E2_cheating)

        MSE_V_array[i]=deepcopy(MSE_V)
        MSE_δ_array[i]=deepcopy(MSE_δ)

        MSE_V_X_array[i]=deepcopy(MSE_V_X)
        MSE_δ_X_array[i]=deepcopy(MSE_δ_X)

        MSE_V_cheating_array[i]=deepcopy(MSE_V_cheating)
        MSE_δ_cheating_array[i]=deepcopy(MSE_δ_cheating)

        MSE_V_X_cheating_array[i]=deepcopy(MSE_V_X_cheating)
        MSE_δ_X_cheating_array[i]=deepcopy(MSE_δ_X_cheating)

        J_arr[i] = deepcopy(J)
        J_X_arr[i] = deepcopy(J_X)

end

#Save generated measurements for any snapshot
snap_shot_number=55
save_measurements(System_array,snap_shot_number)

#Mean Squared Errors
MSE_V=mean(MSE_V_array)
MSE_V_X=mean(MSE_V_X_array)
MSE_δ=mean(MSE_δ_array)
MSE_δ_X=mean(MSE_δ_X_array)

MSE_V_cheating=mean(MSE_V_cheating_array)
MSE_V_X_cheating=mean(MSE_V_X_cheating_array)
MSE_δ_cheating=mean(MSE_δ_cheating_array)
MSE_δ_X_cheating=mean(MSE_δ_X_cheating_array)

geomean(MSE_V)
geomean(MSE_V_X)
geomean(E1_array)
geomean(E2_array)
#Error indices median values
@show median(E1_array)
@show median(E2_array)
@show median(E1_cheating_array)
@show median(E2_cheating_array)

#Increase in accuracy achieved
Accuracy_Increase = (median(E1_array)-median(E2_array))/median(E1_array)

Time = 148
        r, Z,h = calculate_residuals(System_array_Virtual[Time])
        System_array_Virtual[Time].Estimated_wind_Speed
        Ω,Jacobian = calculate_cov_matrix(System_array_Virtual[Time],true)
        B = diag(Ω)
        A = r./sqrt.(B)
        bar(abs.(A),legend = false,xlabel = "Measurement Index",ylabel = "Normalized Residual")
        savefig(string("Normalized Residuals of Time Instance ",string(Time),".pdf"))

r_w = zeros(192,1)

r_w = Vector{Float64}(undef,192)
for Time in 1:192
        r, Z,h = calculate_residuals(System_array_Virtual[Time])
        System_array_Virtual[Time].Estimated_wind_Speed
        try
                Ω,J = calculate_cov_matrix(System_array_Virtual[Time],true)
                B = diag(Ω)
                A = r./sqrt.(B)
                if A[40] >= 100
                        A[40] = 75
                end
                r_w[Time] = A[40]
        catch y
                println(y)
        end
end


r_w[isnan.(r_w) .==0]

r_w[abs.(r_w) .>= norminvcdf(1-0.01/2)]

Limits = abs.(r_w) .>= norminvcdf(1-0.01/2)


TruePositive = zeros(50,1)
TrueNegative = zeros(145,1)
FalsePositive = zeros(15,1)
FalseNegative = zeros(40,1)
Undefined_DRG = zeros(50,1)
Undefined = zeros(50,1)

let a = a ,r_w = r_w

        u_DRG = 1
        u = 1
        TP = 1
        TN = 1
        FN = 1
        FP = 1
        for i in 1:192
                if isnan(r_w[i]) && i <= 144
                        Undefined[u] = i
                        u = u + 1
                elseif isnan(r_w[i]) && i > 144
                        Undefined_DRG[u_DRG] = i
                        u_DRG = u_DRG + 1
                else
                        if Limits[i] == 1 && i> 144
                                TruePositive[TP] = i
                                TP = TP + 1
                        elseif Limits[i] == 1 && i <= 144
                                FalsePositive[FP] = i
                                FP = FP +1
                        elseif Limits[i] == 0 && i>144
                                FalseNegative[FN] = i
                                FN = FN +1
                        elseif Limits[i] == 0 && i<= 144
                                TrueNegative[TN] = i
                                TN = TN +1
                        end
                end

        end
end

TruePositive = TruePositive[TruePositive .!= 0]
TrueNegative = TrueNegative[TrueNegative .!= 0]
FalsePositive = FalsePositive[FalsePositive .!= 0]
FalseNegative = FalseNegative[FalseNegative .!= 0]
Undefined = Undefined[Undefined .!= 0]
Undefined_DRG = Undefined_DRG[Undefined_DRG .!= 0]

bar(abs.(r_w),legend = false,xlabel = "Time (hours)",ylabel = "Normalized Residual", title = string("Normalized Residuals of Exogenous Measurement"))
savefig("Normalized_residuals.pdf")

norminvcdf(1-0.05/2)

#Saving variables for the current run (date: 20 Jan,2020)
df = DataFrame(ErrorIndex_NoCheating_NoX = median(E1_array),ErrorIndex_NoCheating_X = median(E2_array),
        ErrorIndex_Cheating_NoX = median(E1_cheating_array),ErrorIndex_Cheating_X = median(E2_cheating_array),
        Accuracy_Increase_NoCheating = Accuracy_Increase)
CSV.write(string(Dates.format(Dates.now(), "dd-u-yyyy_HH-MM-SS"),".csv"),df)


#Error indices boxplot
StatsPlots.boxplot([E1_array E2_array E1_cheating_array E2_cheating_array],title="Error index boxplot",
        ylabel="Error Index",xticks=(1:4,["E" "E_X" "E_DRG" "E_DRG_X"]),legend=false)
savefig("Error_index_boxplots.pdf")


#Degrees of freedom
k = size(System_array[1].Measurements_PS,1) - (2*System.N_bus-1)
k_x = k+1


#Some informative graphs
ctg = repeat(["w/o exogenous", "with exogenous"], inner = 5)
nam = repeat(1:5, outer = 2)

StatsPlots.groupedbar(nam, [MSE_δ MSE_δ_X], group = ctg, xlabel = "Buses", ylabel = "Mean Squared Error",
        title = "MSE of Phase angles", bar_width = 0.67,
        lw = 0, framestyle = :box)
savefig("MSE_delta.pdf")


StatsPlots.groupedbar(nam, [MSE_V MSE_V_X], group = ctg, xlabel = "Buses", ylabel = "Mean Squared Error",
        title = "MSE of bus Voltages", bar_width = 0.67,
        lw = 0, framestyle = :box)
savefig("MSE_V.pdf")


StatsPlots.groupedbar(nam, [MSE_δ_cheating MSE_δ_X_cheating], group = ctg, xlabel = "Buses", ylabel = "Mean Squared Error",
        title = "MSE of Phase angles in case of cheating", bar_width = 0.67,
        lw = 0, framestyle = :box)
savefig("MSE_delta_cheating.pdf")

StatsPlots.groupedbar(nam, [MSE_V_cheating MSE_V_X_cheating], group = ctg, xlabel = "Buses", ylabel = "Mean Squared Error",
        title = "MSE of bus Voltages in case of cheating", bar_width = 0.67,
        lw = 0, framestyle = :box)
savefig("MSE_V_cheating.pdf")


plot([chisqccdf(k,median([J_arr[j][i][1] for j in 1:n])) for i in 1:192])
plot!([chisqccdf(k,median([J_X_arr[j][i][1] for j in 1:n])) for i in 1:192])
savefig("Confidence_Curve_max_min_med.pdf")
