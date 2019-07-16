mutable struct p_base
    k::Array{Float64,1};
    x_o2aq::Float64;
    x_h2o::Float64;
    x_h2o2::Float64;
end

mutable struct p_parallel
    #All old stuff
    k::Array{Float64,1};
    x_o2aq::Float64;
    x_h2o::Float64;
    x_h2o2::Float64;

    #Plus some new stuff
    U::Float64#Voltage
    delta_G::Array{Float64,1}#Delta G's
    chem_E::Array{Float64,1} #E's for chemical steps
    A::Float64#A
    E::Float64
    beta::Float64
    kb::Float64
    T::Float64
    h::Float64
end
