using DifferentialEquations
using Flux
using LinearAlgebra
using Random
using Distributions

include("kinetic_rate_out.jl")
include("calc_ks.jl")
include("generate_test_data.jl")
include("rate_equations.jl")

num_points=500


#The initial input vector will be the dict of energies
model = Chain(Dense(13,50),Dense(50,50),Dense(50,10),softmax)
loss(x,y)=Flux.mse(model(x),y)

x_data,y_data = generate_test_data(num_points)
data=[(x_data,y_data)]
ps=Flux.params(model)

cb() = println("Training!")

Flux.train!(loss,ps,data,ADAM(0.001),cb = () -> println("training"))
