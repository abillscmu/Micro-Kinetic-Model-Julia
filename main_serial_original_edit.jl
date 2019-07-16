using DifferentialEquations
using Sundials
using Plots
include("rate_equations.jl")
include("datastructures.jl")
function main_serial(;bigPlot=false,smallPlot=false)


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
    h=4.14e-15
    #Calculation of Rate Constants
    k_pos = zeros(13)
    K = zeros(13)

    #Constant Constants
    kb = 8.617 .* (10 .^ -5)
    T = 298
    kbT = kb .* T
    kbTh = kbT ./ h

    #Trib
    E = 0.26 + 0.05
    beta = 0.5

    delta_G_1 = 0
    delta_G_2 = -0.2
    delta_G_3 = -.14
    delta_G_4 = -1.35
    delta_G_5 = -0.05
    delta_G_6 = 0.15
    delta_G_7 = -.93
    delta_G_8 = -.38
    delta_G_9 = -0.2
    delta_G_10 = -0.83
    delta_G_11 = 0.28
    delta_G_12 = -1.49
    delta_G_13 = 0.32


    G_0_3 = -0.9+delta_G_3
    G_0_4 = -0.9+delta_G_4
    G_0_5 = -0.9+delta_G_5
    G_0_6 = -0.9+delta_G_6
    G_0_8 = -0.9+delta_G_8
    G_0_9 = -0.9+delta_G_9
    G_0_11 = -0.9+delta_G_11

    E_1 = 0
    E_2 = 0
    E_7 = 0.48
    E_10 = 0.37
    E_12 = 0.46
    E_13 = 0

    #Delta G error 0.05
    #lognormal exponential of Gaussian

    U_3 = -G_0_3
    U_4 = -G_0_4
    U_5 = -G_0_5
    U_6 = -G_0_6
    U_8 = -G_0_8
    U_9 = -G_0_9
    U_11 = -G_0_11
    A = 1.0 .* (10 .^ 9)



    #Potential: Input
for n = 1:num_points
        U = U_vec[n]

        k_pos[1] = (8 .* (10 .^ 5)) .* min(1,exp.(-E_1 ./ kbT ))
        k_pos[2] = (1 .* (10 .^ 8)) .* exp.(-E_2 ./ kbT )
        k_pos[3] = A .* min(1,exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_3))) ./ kbT))
        k_pos[4] = A .* min(1,exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_4))) ./ kbT))
        k_pos[5] = A .*min(1, exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_5))) ./ kbT))
        k_pos[6] = A .* min(1,exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_6))) ./ kbT))
        k_pos[7] = kbTh .* min(1,exp.(-E_7 ./ kbT))
        k_pos[8] = A .* min(1,exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_8))) ./ kbT))
        k_pos[9] = A .* min(1,exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_9))) ./ kbT))
        k_pos[10] = kbTh .* min(1,exp.(-E_10 ./ kbT))
        k_pos[11] = A .* min(1,exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_11))) ./ kbT))
        k_pos[12] = kbTh .* min(1,exp.(-E_12./ kbT))
        k_pos[13] = (1.00 .* (10 .^ 8)) .* min(1,exp.(-E_13 ./ kbT))

        K[1] = exp.(- (delta_G_1 ./ kbT))
        K[2] = exp.(- (delta_G_2 ./ kbT))
        K[3] = exp.(- ((G_0_3+U) ./ kbT))
        K[4] = exp.(- ((G_0_4+U) ./ kbT))
        K[5] = exp.(- ((G_0_5+U) ./ kbT))
        K[6] = exp.(- ((G_0_6+U) ./ kbT))
        K[7] = exp.(- (delta_G_3 ./ kbT))
        K[8] = exp.(- ((G_0_8+U) ./ kbT))
        K[9] = exp.(- ((G_0_9+U) ./ kbT))
        K[10] = exp.(- (delta_G_10 ./ kbT))
        K[11] = exp.(- ((G_0_11+U) ./ kbT))
        K[12] = exp.(- (delta_G_12 ./ kbT))
        K[13] = exp.(- (delta_G_13 ./ kbT))

        k_neg = k_pos ./ K

        k_init = vcat(k_pos,k_neg)
        #k_init = ones(26)
        x_o2aq = 2.34 .* (10 .^ -5)
        x_h2o = 1
        x_h2o2 = 0
        #for n = 1:26
        #        if abs(k_init[n]>1e20)
        #            k_init[n]=1e20
        #        end
        #    end

        #Creation of parameter object
        p = p_base(k_init,x_o2aq,x_h2o,x_h2o2);

        #Problem setup & Initial Conditions
        diff_variables = trues(10);
        diff_variables[7]=false;
        diff_variables[10]=false

        tspan = (0,2000.)
        #println(k_init[11])
        #println(k_init[24])

        u_0 = [0.01,1.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]
        du_0 = [0.0,.1,.1,.1,.1,.1,0.1,.1,0.1,0.1]

        prob = ODEProblem(rate_equations,u_0,tspan,p)
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

end
