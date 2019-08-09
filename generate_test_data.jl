using Random
include("calc_ks.jl")
function generate_test_data(num_points)
test_data_x = Array{Float32,2}(undef,13,num_points)
test_data_y = Array{Float32,2}(undef,10,num_points)

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
G_U0s = Dict{String,Float32}("O2aq"=>1.32+dGO2solv,"O2dl"=>1.32+dGO2solv,"O2"=>1.392+(1.2*dGOH),"OOH"=>1.25+dGOH,"O"=>-0.14+0.04+2*dGOH,"OH"=>-0.15+dGOH,"HOOH"=>1.533+0.04*dGOH,"Ob"=>0.5670+(2*dGOH),"OHb"=>0.1979+dGOH,"H2O2aq"=>1.75+H2O2aq_corr ," "=>0.,"b"=>0.)
qs = Dict{String,Float32}("O2aq"=>-4.,"O2dl"=>-4.,"O2"=>-4.,"OOH"=>-3.,"O"=>-2.,"OH"=>-1.,"HOOH"=>-2.,"Ob"=>-2.,"OHb"=>-1.,"H2O2aq"=>-2.," "=>0.,"b"=>0.)
G_Us = Dict{String,Float32}("O2aq"=>1.32+dGO2solv,"O2dl"=>1.32+dGO2solv,"O2"=>1.392+(1.2*dGOH),"OOH"=>1.25+dGOH,"O"=>-0.14+0.04+2*dGOH,"OH"=>-0.15+dGOH,"HOOH"=>1.533+0.04*dGOH,"Ob"=>0.5670+(2*dGOH),"OHb"=>0.1979+dGOH,"H2O2aq"=>1.75+H2O2aq_corr ," "=>0,"b"=>0.)
rng=MersenneTwister(1234)
for state in keys(G_U0s)
    G_U0s[state] = G_U0s[state]+(qs[state]*(-0.9));
end

g_0s = [G_U0s["O2aq"],G_U0s["O2dl"],G_U0s["O2"],G_U0s["OOH"],G_U0s["O"],G_U0s["OH"],G_U0s["HOOH"],G_U0s["Ob"],G_U0s["OHb"],G_U0s["H2O2aq"],G_U0s[" "],G_U0s["b"]]
    #Potential: Input
for n=1:num_points
    voltage = 0.2f0+0.8f0*rand(Float32)
    g_vec = vcat(voltage,g_0s[1:10] .+ .05f0*randn(rng,Float32,10),g_0s[11:12])
    test_data_x[:,n]=g_vec
    test_data_y[:,n]=calc_ans(g_vec)
end
return test_data_x,test_data_y
end
