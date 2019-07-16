using DifferentialEquations
using JLD2

@load "test.jld2" sim_arr
num_points = length(sim_arr)

high=Array{Any,1}(undef,num_points)
low=Array{Any,1}(undef,num_points)



for n = 1:num_points
    sim_this_round = sim_arr[n]
    high[n] = sim_this_round.qhigh[end]
    low[n] = sim_this_round.qlow[end]

end
