using DifferentialEquations
using Random
using Sundials
using DelimitedFiles
using Distributions
include("rate_equations.jl")
include("datastructures.jl")
include("kinetic_rate_mc_run.jl")
using SharedArrays

num_runs=1000
num_points=50
o2dl=Array{Float64,2}(undef,num_points,num_runs)
stara=Array{Float64,2}(undef,num_points,num_runs)
o2stara=Array{Float64,2}(undef,num_points,num_runs)
theta_oohstarA=Array{Float64,2}(undef,num_points,num_runs)
theta_ostarA=Array{Float64,2}(undef,num_points,num_runs)
theta_ohstarA=Array{Float64,2}(undef,num_points,num_runs)
theta_h2o2starA=Array{Float64,2}(undef,num_points,num_runs)
theta_ohstarB=Array{Float64,2}(undef,num_points,num_runs)
theta_starB=Array{Float64,2}(undef,num_points,num_runs)
theta_ostarB=Array{Float64,2}(undef,num_points,num_runs)
timetrack=Array{Float64,2}(undef,num_points,num_runs)



for r = 1:num_runs
	o2dl[:,r],stara[:,r],o2stara[:,r],theta_oohstarA[:,r],theta_ostarA[:,r],theta_ohstarA[:,r],theta_h2o2starA[:,r],theta_ohstarB[:,r],theta_starB[:,r],theta_ostarB[:,r],timetrack[:,r]=kinetic_rate_mc_run(r)
	println(r)
end
