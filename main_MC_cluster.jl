@everywhere using DifferentialEquations
@everywhere using Random
@everywhere using Sundials
@everywhere using LinearAlgebra
using DelimitedFiles
@everywhere include("rate_equations.jl")
@everywhere include("calc_ks.jl")
using SharedArrays

@everywhere num_runs=1000
@everywhere num_points=50

@everywhere include("initialize.jl")

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
timetrack=SharedArray{Float64,2}(num_points,num_runs)

#Potential: Input
@sync @distributed for r=1:num_runs
for n = 1:num_points
		start_time = time_ns()
    	U0=0.9
		U=U_vec[n]
		k_init=calc_ks(vcat(U,g_U0s),mc=true,seed=r)
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
        p = vcat(k_init,x_o2aq,x_h2o,x_h2o2);

        #Problem setup & Initial Conditions
        mm=ones(10)
		mm[7]=0
		mm[10]=0
		mm=Diagonal(mm)

        tspan = (0,0.1)
        #println(k_init[11])
        #println(k_init[24])

        u_0 = [0.01,1.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]
        du_0 = [0.0,.1,.1,.1,.1,.1,0.1,.1,0.1,0.1]

		func = ODEFunction(rate_equations,mass_matrix=mm)
        prob = ODEProblem(func,u_0,tspan,p)
        sol = solve(prob,Rodas4())

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

	    timetrack[n,r]=(time_ns()-start_time)/1e9




    end
        if(r%100==0)
            println("vvvITERATIONvvv")
            println(r)
        end

end



open("newdata/o2dl.csv","w") do io
    writedlm(io,o2dl,',');
end

open("newdata/stara.csv","w") do io
    writedlm(io,stara,',')
end

open("newdata/o2stara.csv","w") do io
    writedlm(io,o2stara,',')
end

open("newdata/theta_oohstarA.csv","w") do io
    writedlm(io,theta_oohstarA,',')
end

open("newdata/theta_ostarA.csv","w") do io
    writedlm(io,theta_ostarA,',')
end

open("newdata/theta_ohstarA.csv","w") do io
    writedlm(io,theta_ohstarA,',')
end

open("newdata/theta_h2o2starA.csv","w") do io
    writedlm(io,theta_h2o2starA,',')
end

open("newdata/theta_ohstarB.csv","w") do io
    writedlm(io,theta_ohstarB,',')
end

open("newdata/theta_starB.csv","w") do io
    writedlm(io,theta_starB,',')
end

open("newdata/theta_ostarB.csv","w") do io
    writedlm(io,theta_ostarB,',')
end

open("newdata/timetrack.csv","w") do io
   writedlm(io,timetrack,',')
end
