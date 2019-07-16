using Distributed
using JLD2

@everywhere include("rate_equations.jl")
@everywhere include("datastructures.jl")
@everywhere include("mc_modifier_fake.jl")
@everywhere using DifferentialEquations
@everywhere using Sundials

function main_MC(;bigPlot=false,smallPlot=false)


num_points=10
sim_arr=Array{Any,1}(undef,num_points)
delta_G = zeros(13)
chem_E = zeros(13)

U_vec = range(0.5,stop=1,length=num_points);

#Constant Constants
h=4.14e-15
kb = 8.617 .* (10 .^ -5)
T = 298
kbT = kb .* T
kbTh = kbT ./ h


#Trib
E = 0.26
beta = 0.5

delta_G[1] = 0
delta_G[2]  = -0.2
delta_G[3]  = -.14
delta_G[4] = -1.35
delta_G[5]  = -0.05
delta_G[6]  = 0.15
delta_G[7]  = -.93
delta_G[8]  = -.38
delta_G[9]  = -0.2
delta_G[10]  = -0.83
delta_G[11]  = 0.28
delta_G[12]  = -1.49
delta_G[13]  = 0.32



chem_E[1] = 0
chem_E[2] = 0
chem_E[3] = 0.48
chem_E[4] = 0.37
chem_E[5] = 0.46
chem_E[6] = 0

A = 1.0 .* (10 .^ 9)

    #Potential: Input
@time for n = 1:num_points
        U = U_vec[n]


        k_init = zeros(26)
        #k_init = ones(26)
        x_o2aq = 2.34 .* (10 .^ -5)
        x_h2o = 1
        x_h2o2 = 0

        p = p_parallel(k_init,x_o2aq,x_h2o,x_h2o2,U,delta_G,chem_E,A,E,beta,kb,T,h);

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
        monte_prob = MonteCarloProblem(prob,prob_func=mc_modifier)
        sim = solve(monte_prob,IDA(max_num_iters_ic = 1000),MonteDistributed(),num_monte=10000)
        sim_arr[n] = MonteCarloSummary(sim)

        println("Done with Iteration")
        println(n)

end

@save "test_fake.jld2" sim_arr
end
