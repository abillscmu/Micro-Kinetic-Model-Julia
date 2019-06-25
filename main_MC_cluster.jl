@everywhere using DifferentialEquations
@everywhere using Sundials
using DelimitedFiles
@everywhere include("rate_equations.jl")
@everywhere include("datastructures.jl")
using SharedArrays

@everywhere num_runs = 10000
@everywhere num_points=50
o2dl=SharedArray{Float64,2}(num_points,num_runs)
stara=SharedArray{Float64,2}(num_points,num_runs)
o2stara=SharedArray{Float64,2}(num_points,num_runs)
theta_oohstarA=SharedArray{Float64,2}(num_points,num_runs)
theta_ostarA=SharedArray{Float64,2}(num_points,num_runs)
theta_ohstarA=SharedArray{Float64,2}(num_points,num_runs)
theta_h2o2starA=SharedArray{Float64,2}(num_points,num_runs)
theta_ohstarB=SharedArray{Float64,2}(num_points,num_runs)
theta_starB=SharedArray{Float64,2}(num_points,num_runs)
theta_ostarB=SharedArray{Float64,2}(num_points,num_runs)

@sync @distributed for r = 1:num_runs




    #Potential: Input
        for n = 1:num_points
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
        delta_G_2 = -0.2 + 0.05
        delta_G_3 = -.14 + 0.05
        delta_G_4 = -1.35 + 0.05
        delta_G_5 = -0.05 + 0.05
        delta_G_6 = 0.15 + 0.05
        delta_G_7 = -.93 + 0.05
        delta_G_8 = -.38 + 0.05
        delta_G_9 = -0.2 + 0.05
        delta_G_10 = -0.83 + 0.05
        delta_G_11 = 0.28 + 0.05
        delta_G_12 = -1.49 + 0.05
        delta_G_13 = 0.32 + 0.05


        G_0_3 = -0.9+delta_G_3
        G_0_4 = -0.9+delta_G_4
        G_0_5 = -0.9+delta_G_5
        G_0_6 = -0.9+delta_G_6
        G_0_8 = -0.9+delta_G_8
        G_0_9 = -0.9+delta_G_9
        G_0_11 = -0.9+delta_G_11

        E_1 = 0
        E_2 = 0
        E_7 = 0.48 + 0.05
        E_10 = 0.37 + 0.05
        E_12 = 0.46 + 0.05
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

        U = U_vec[n]

        k_pos[1] = (8 .* (10 .^ 5)) .* exp.(-E_1 ./ kbT )
        k_pos[2] = (1 .* (10 .^ 8)) .* exp.(-E_2 ./ kbT )
        k_pos[3] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_3))) ./ kbT)
        k_pos[4] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_4))) ./ kbT)
        k_pos[5] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_5))) ./ kbT)
        k_pos[6] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_6))) ./ kbT)
        k_pos[7] = kbTh .* exp.(-E_7 ./ kbT)
        k_pos[8] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_8))) ./ kbT)
        k_pos[9] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_9))) ./ kbT)
        k_pos[10] = kbTh .* exp.(-E_10 ./ kbT)
        k_pos[11] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_11))) ./ kbT)
        k_pos[12] = kbTh .* exp.(-E_12./ kbT)
        k_pos[13] = (1.00 .* (10 .^ 8)) .* exp.(-E_13 ./ kbT)

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

        tspan = (0,0.1)
        #println(k_init[11])
        #println(k_init[24])

        u_0 = [0.01,1.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]
        du_0 = [0.0,.1,.1,.1,.1,.1,0.1,.1,0.1,0.1]

        prob = DAEProblem(rate_equations,du_0,u_0,tspan,p,differential_vars=diff_variables)
        sol = solve(prob,IDA(max_num_iters_ic = 10000))

        o2dl[n,r]=sol.u[end][1]
        stara[n,r]=sol.u[end][2]
        o2stara[n,r]=sol.u[end][7]
        theta_oohstarA[n,r]=sol.u[end][4]
        theta_ostarA[n,r]=sol.u[end][5]
        theta_ohstarA[n,r]=sol.u[end][6]
        theta_h2o2starA[n,r]=sol.u[end][3]

        theta_ohstarB[n,r]=sol.u[end][8]
        theta_starB[n,r]=sol.u[end][10]
        theta_ostarB[n,r]=sol.u[end][9]

    end
        if(r%100==0)
            println("vvvITERATIONvvv")
            println(r)
        end

end


open("o2dl.csv","w") do io
    writedlm(io,o2dl,',');
end

open("stara.csv","w") do io
    writedlm(io,stara,',')
end

open("o2stara.csv","w") do io
    writedlm(io,o2stara,',')
end

open("theta_oohstarA.csv","w") do io
    writedlm(io,theta_oohstarA,',')
end

open("theta_ostarA.csv","w") do io
    writedlm(io,theta_ostarA,',')
end

open("theta_ohstarA.csv","w") do io
    writedlm(io,theta_ohstarA,',')
end

open("theta_h2o2starA.csv","w") do io
    writedlm(io,theta_h2o2starA,',')
end

open("theta_ohstarB.csv","w") do io
    writedlm(io,theta_ohstarB,',')
end

open("theta_starB.csv","w") do io
    writedlm(io,theta_starB,',')
end

open("theta_ostarB.csv","w") do io
    writedlm(io,theta_ostarB,',')
end
