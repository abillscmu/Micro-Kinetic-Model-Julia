using DifferentialEquations
using Flux
using LinearAlgebra
using Random
using Distributions
using Statistics

include("kinetic_rate_out.jl")
include("calc_ks.jl")
include("generate_test_data.jl")
include("rate_equations.jl")

#Number of data points on which to train
num_points=100000


#Create Model
model = Chain(Dense(13,12,tanh),Dense(12,12,tanh),Dense(12,11,tanh),Dense(11,10,tanh))

#Define Loss Function
loss(x,y)=Flux.mse(model(x),y)

#Generate Test Data
x_data,y_data = generate_test_data(num_points)

#Condition Data for NN

#Log y data
y_data=log10.(y_data)

#Standardize Outputs
y_data_standardized = similar(y_data)

#Preallocate to save later
y_mean_vec = zeros(Float32,size(y_data)[1])
y_std_vec = ones(Float32,size(y_data)[1])

for n = 1:size(y_data)[1]
    y_mean_vec[n] = mean(y_data[n,:])
    y_std_vec[n] = std(y_data[n,:])
    y_data_standardized[n,:] = (y_data[n,:] .- y_mean_vec[n]) ./ y_std_vec[n];
end

#Standardize Inputs
x_data_standardized = similar(x_data)

#need to preallocate to save later
x_mean_vec=zeros(Float32,size(x_data)[1])
x_std_vec = ones(Float32,size(x_data)[1])

for n = 1:size(x_data)[1]
    if(n<12)
        x_mean_vec[n]=mean(x_data[n,:])
        x_std_vec[n] = std(x_data[n,:])
    end
    x_data_standardized[n,:] = (x_data[n,:] .- x_mean_vec[n]) ./ x_std_vec[n]
end



data=[(x_data_standardized,y_data_standardized)]
ps=Flux.params(model)

evalcb() = @show(loss(x_data_standardized[:,25],y_data_standardized[:,25]))

using Flux: @epochs, throttle
@epochs 500 Flux.train!(loss,ps,data,ADAM(0.001),cb = throttle(evalcb,5))



using Plots

include("rate_equations.jl")
include("calc_ks.jl")
bigPlot=true
smallPlot=true


num_points=50
o2dl=Array{Float64,1}(undef,num_points)
stara=Array{Float64,1}(undef,num_points)
o2stara=Array{Float64,1}(undef,num_points)
theta_oohstarA=Array{Float64,1}(undef,num_points)
theta_ostarA=Array{Float64,1}(undef,num_points)
theta_ohstarA=Array{Float64,1}(undef,num_points)
theta_h2o2starA=Array{Float64,1}(undef,num_points)
theta_ohstarB=Array{Float64,1}(undef,num_points)
theta_starB=Array{Float64,1}(undef,num_points)
theta_ostarB=Array{Float64,1}(undef,num_points)

o2dl_surr=Array{Any,1}(undef,num_points)
stara_surr=Array{Any,1}(undef,num_points)
o2stara_surr=Array{Any,1}(undef,num_points)
theta_oohstarA_surr=Array{Any,1}(undef,num_points)
theta_ostarA_surr=Array{Any,1}(undef,num_points)
theta_ohstarA_surr=Array{Any,1}(undef,num_points)
theta_h2o2starA_surr=Array{Any,1}(undef,num_points)
theta_ohstarB_surr=Array{Any,1}(undef,num_points)
theta_starB_surr=Array{Any,1}(undef,num_points)
theta_ostarB_surr=Array{Any,1}(undef,num_points)


U_vec = range(0.2,stop=1,length=num_points);
    #Calculation of Rate Constants
k_pos = zeros(13)
K = zeros(13)


    #Constant Constants
kb = 8.617 .* (10 .^ -5)
T = 298
h=4.14e-15
beta = 0.5
kbT = kb .* T
kbTh = kbT ./ h
E = 0.26

cH2O = 1e3/(16. + (2*1.008))
kH_H2O = 1.3e-3
dGO2solv = kbT*log(cH2O/kH_H2O)
dGH2O2solv = kbT*log(cH2O)
dGOH=0
U0 = 0.9 #Potential at which other energies are defined
OHB_destabilization = 0.26+0.1 #*b destabilization energies (Tripkovic et al)
Ob_destabilization = 0.0#see above
H2O2aq_corr = dGH2O2solv
O2ads_corr=0.0
G_U0s = Dict("O2aq"=>1.32+dGO2solv,"O2dl"=>1.32+dGO2solv,"O2"=>1.392+(1.2*dGOH),"OOH"=>1.25+dGOH,"O"=>-0.14+0.04+2*dGOH,"OH"=>-0.15+dGOH,"HOOH"=>1.533+0.04*dGOH,"Ob"=>0.5670+(2*dGOH),"OHb"=>0.1979+dGOH,"H2O2aq"=>1.75+H2O2aq_corr ," "=>0.,"b"=>0.)
qs = Dict("O2aq"=>-4.,"O2dl"=>-4.,"O2"=>-4.,"OOH"=>-3.,"O"=>-2.,"OH"=>-1.,"HOOH"=>-2.,"Ob"=>-2.,"OHb"=>-1.,"H2O2aq"=>-2.," "=>0.,"b"=>0.)
G_Us = Dict("O2aq"=>1.32+dGO2solv,"O2dl"=>1.32+dGO2solv,"O2"=>1.392+(1.2*dGOH),"OOH"=>1.25+dGOH,"O"=>-0.14+0.04+2*dGOH,"OH"=>-0.15+dGOH,"HOOH"=>1.533+0.04*dGOH,"Ob"=>0.5670+(2*dGOH),"OHb"=>0.1979+dGOH,"H2O2aq"=>1.75+H2O2aq_corr ," "=>0,"b"=>0.)

for state in keys(G_U0s)
    G_U0s[state] = G_U0s[state]+(qs[state]*(-0.9));
end

g_0s = [G_U0s["O2aq"],G_U0s["O2dl"],G_U0s["O2"],G_U0s["OOH"],G_U0s["O"],G_U0s["OH"],G_U0s["HOOH"],G_U0s["Ob"],G_U0s["OHb"],G_U0s["H2O2aq"],G_U0s[" "],G_U0s["b"]]
    #Potential: Input
for n = 1:num_points
        U = U_vec[n]


        ground_truth=calc_ans(vcat(U,g_0s))


        o2dl[n]=ground_truth[1]
        stara[n]=ground_truth[2]
        o2stara[n]=ground_truth[7]
        theta_oohstarA[n]=ground_truth[4]
        theta_ostarA[n]=ground_truth[5]
        theta_ohstarA[n]=ground_truth[6]
        theta_h2o2starA[n]=ground_truth[3]

        theta_ohstarB[n]=ground_truth[8]
        theta_starB[n]=ground_truth[10]
        theta_ostarB[n]=ground_truth[9]

        test_x_data = vcat(U,g_0s)
        test_x_data_standardized = (test_x_data .- x_mean_vec) ./ x_std_vec
        nn_prediction=Tracker.data(model(test_x_data_standardized))
        #Destandardize Y Data
        nn_prediction = (nn_prediction .* y_std_vec) .+ y_mean_vec

        o2dl_surr[n]=10 .^ nn_prediction[1]
        stara_surr[n]=10 .^ nn_prediction[2]
        o2stara_surr[n]=10 .^ nn_prediction[7]
        theta_oohstarA_surr[n]=10 .^ nn_prediction[4]
        theta_ostarA_surr[n]=10 .^ nn_prediction[5]
        theta_ohstarA_surr[n]=10 .^ nn_prediction[6]
        theta_h2o2starA_surr[n]=10 .^ nn_prediction[3]

        theta_ohstarB_surr[n]=10 .^ nn_prediction[8]
        theta_starB_surr[n]=10 .^ nn_prediction[10]
        theta_ostarB_surr[n]=10 .^ nn_prediction[9]

end

if(bigPlot)
    plot(U_vec,o2dl,yscale=:log10,lw=3,label="o2dl",color=:green)
    plot!(U_vec,stara,yscale=:log10,lw=3,label="*a",color=:red)
    plot!(U_vec,o2stara,yscale=:log10,lw=3,label="o2*a",color=:blue)
    plot!(U_vec,theta_oohstarA,yscale=:log10,lw=3,label="ooh*a",color=:cyan)
    plot!(U_vec,theta_ostarA,yscale=:log10,lw=3,label="o*a",color=:purple)
    plot!(U_vec,theta_ohstarA,yscale=:log10,lw=3,label="oh*a",color=:gold)
    plot!(U_vec,theta_h2o2starA,yscale=:log10,lw=3,label="h2o2*a",legend=:bottomleft,color=:black)

    plot!(U_vec,o2dl_surr,yscale=:log10,lw=3,label="o2dl",color=:green,linestyle = :dot)
    plot!(U_vec,stara_surr,yscale=:log10,lw=3,label="*a",color=:red,linestyle = :dot)
    plot!(U_vec,o2stara_surr,yscale=:log10,lw=3,label="o2*a",color=:blue,linestyle = :dot)
    plot!(U_vec,theta_oohstarA_surr,yscale=:log10,lw=3,label="ooh*a",color=:cyan,linestyle = :dot)
    plot!(U_vec,theta_ostarA_surr,yscale=:log10,lw=3,label="o*a",color=:purple,linestyle = :dot)
    plot!(U_vec,theta_ohstarA_surr,yscale=:log10,lw=3,label="oh*a",color=:gold,linestyle = :dot)
    plot!(U_vec,theta_h2o2starA_surr,yscale=:log10,lw=3,label="h2o2*a",legend=:bottomleft,color=:black,linestyle = :dot)
    savefig("BigPlot.png")
end
if(smallPlot)
    plot(U_vec,stara,lw=3,label="*a",color=:red)
    plot!(U_vec,theta_ostarA,lw=3,label="o*a",color=:purple)
    plot!(U_vec,theta_ohstarA,lw=3,label="oh*a",color=:gold)

    plot!(U_vec,stara_surr,lw=3,label="*a",color=:red,linestyle = :dot)
    plot!(U_vec,theta_ostarA_surr,lw=3,label="o*a",color=:purple,linestyle = :dot)
    plot!(U_vec,theta_ohstarA_surr,lw=3,label="oh*a",color=:gold,linestyle = :dot)
    savefig("SmallPlot.png")
end
#Parameter Structure: [G_U0s, Voltage]

#The initial input vector will be the dict of energies
