using DifferentialEquations
using Sundials
using Plots
using LinearAlgebra
include("rate_equations.jl")
include("calc_ks.jl")
bigPlot=true
smallPlot=false


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


        k_init=calc_ks(vcat(U,g_0s))
        #k_init = ones(26)
        x_o2aq = 2.34 .* (10 .^ -5)
        x_h2o = 1
        x_h2o2 = 0


        #Creation of parameter object
        p = vcat(k_init,x_o2aq,x_h2o,x_h2o2);

        #Problem setup & Initial Conditions
        diff_variables = trues(10);
        diff_variables[7]=false;
        diff_variables[10]=false

        tspan = (0,10.)

        u_0 = [1.,0.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]
        du_0 = [0.0,.1,-.1,.1,.1,.1,0.1,.1,0.1,0.1]
        mm = ones(10);
        mm[7]=0
        mm[10]=0
        MM=Diagonal(mm);
        func = ODEFunction(rate_equations,mass_matrix=MM)
        prob = ODEProblem(func,u_0,tspan,p)
        sol = solve(prob,Rodas4())

        o2dl[n]=sol.u[end][1]
        stara[n]=sol.u[end][2]
        o2stara[n]=sol.u[end][7]
        theta_oohstarA[n]=sol.u[end][4]
        theta_ostarA[n]=sol.u[end][5]
        theta_ohstarA[n]=sol.u[end][6]
        theta_h2o2starA[n]=sol.u[end][3]

        theta_ohstarB[n]=sol.u[end][8]
        theta_starB[n]=sol.u[end][10]
        theta_ostarB[n]=sol.u[end][9]

end

if(bigPlot)
    plot(U_vec,o2dl,yscale=:log10,lw=3,label="o2dl",color=:green)
    plot!(U_vec,stara,yscale=:log10,lw=3,label="*a",color=:red)
    plot!(U_vec,o2stara,yscale=:log10,lw=3,label="o2*a",color=:blue)
    plot!(U_vec,theta_oohstarA,yscale=:log10,lw=3,label="ooh*a",color=:cyan)
    plot!(U_vec,theta_ostarA,yscale=:log10,lw=3,label="o*a",color=:purple)
    plot!(U_vec,theta_ohstarA,yscale=:log10,lw=3,label="oh*a",color=:gold)
    plot!(U_vec,theta_h2o2starA,yscale=:log10,lw=3,label="h2o2*a",legend=:bottomleft,color=:black)
    savefig("BigPlot.png")
end
if(smallPlot)
    plot(U_vec,stara,lw=3,label="*a",color=:red)
    plot!(U_vec,theta_ostarA,lw=3,label="o*a",color=:purple)
    plot!(U_vec,theta_ohstarA,lw=3,label="oh*a",color=:gold)
    savefig("SmallPlot.png")
end
#Parameter Structure: [G_U0s, Voltage]

#The initial input vector will be the dict of energies
model = Chain(x -> x.^2,
             Dense(16,50,tanh),
             Dense(50,50,tanh),
             Dense(50,10))
