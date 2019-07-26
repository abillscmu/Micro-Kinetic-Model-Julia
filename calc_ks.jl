
using Random
function calc_ks(params3;mc=false,seed=0)
U = params3[1]
if(mc)
    rng=MersenneTwister(seed)
    rnd=0.05 .* randn(rng,12)
    params3[2:end]=params3[2:end] .+ rnd;
end



indices = Dict("O2aq"=>1,"O2dl"=>2,"O2"=>3,"OOH"=>4,"O"=>5,"OH"=>6,"HOOH"=>7,"Ob"=>8,"OHb"=>9,"H2O2aq"=>10 ," "=>11,"b"=>12)
G_U0s  = Dict("O2aq"=>1.32+dGO2solv,"O2dl"=>1.32+dGO2solv,"O2"=>1.392+(1.2*dGOH),"OOH"=>1.25+dGOH,"O"=>-0.14+0.04+2*dGOH,"OH"=>-0.15+dGOH,"HOOH"=>1.533+0.04*dGOH,"Ob"=>0.5670+(2*dGOH),"OHb"=>0.1979+dGOH,"H2O2aq"=>1.75+H2O2aq_corr ," "=>0.,"b"=>0.)
qs = Dict("O2aq"=>-4.,"O2dl"=>-4.,"O2"=>-4.,"OOH"=>-3.,"O"=>-2.,"OH"=>-1.,"HOOH"=>-2.,"Ob"=>-2.,"OHb"=>-1.,"H2O2aq"=>-2.," "=>0.,"b"=>0.)
###ORDERING
G_Us=copy(G_U0s)

for state in keys(indices)
    G_U0s[state]=params3[indices[state]+1]
    G_Us[state] = G_U0s[state]+(qs[state]*U);
end
delta_G_1 = G_Us["O2dl"]-G_Us["O2aq"]
delta_G_2 = G_Us["O2"]-G_Us["O2dl"]
delta_G_3 = G_Us["OOH"]-G_Us["O2"]
delta_G_4 = G_Us["O"]-G_Us["OOH"]
delta_G_5 = G_Us["OH"]-G_Us["O"]
delta_G_6 = G_Us[" "]-G_Us["OH"]
delta_G_7 = G_Us["O"]+G_Us["Ob"]-G_Us["O2"]
delta_G_8 = G_Us["OHb"]-G_Us["Ob"]
delta_G_9 = G_Us["b"]-G_Us["OHb"]
delta_G_10 = G_Us["Ob"]+G_Us["OH"]-G_Us["OOH"]
delta_G_11 = G_Us["HOOH"]-G_Us["OOH"]
delta_G_12 = G_Us["OH"]+G_Us["OHb"]-G_Us["HOOH"]
delta_G_13 = G_Us["H2O2aq"]-G_Us["HOOH"]

beta = 1. / ( kbT )

f0 = exp.(-beta * 0.2244)
f1 = 1
s2 = exp.(-beta * 0.05916)
s3 = 1.
s4 = 1.
s13 = s2

k_pos[1] = (8 .* (10 .^ 5))
k_pos[2] = kbTh *f0 *s2 * min(1.,exp(-beta*delta_G_2))
k_pos[3] = kbTh * f0 * s3 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_3 ) ))
k_pos[4] = kbTh * f0 * s4 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_4 ) ))
k_pos[5] = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_5 ) ))
k_pos[6] = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_6 ) ))
k_pos[7] = kbTh * min(1.,exp( - beta * ( 0.48 + 0.69 * ( 0.9255 + delta_G_7) ) ))
k_pos[8] = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_8 ) ))
k_pos[9] = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_9 ) ))
k_pos[10] = kbTh * min(1.,exp( - beta * ( 0.37 + 0.39 * (0.8330 + delta_G_10) ) ))
k_pos[11] = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_11 ) ))
k_pos[12] =  kbTh * min(1.,exp( - beta * ( 0.462 + 0.19 * (1.4853 + delta_G_12) ) ))
k_pos[13] = kbTh * s13 * f0 * min(1.,exp( - beta * ( delta_G_13 ) ))

K[1] = exp.(- (delta_G_1 ./ kbT))
K[2] = exp.(- (delta_G_2 ./ kbT))
K[3] = exp.(- ((delta_G_3) ./ kbT))
K[4] = exp.(- ((delta_G_4) ./ kbT))
K[5] = exp.(- ((delta_G_5) ./ kbT))
K[6] = exp.(- ((delta_G_6) ./ kbT))
K[7] = exp.(- (delta_G_7 ./ kbT))
K[8] = exp.(- ((delta_G_8) ./ kbT))
K[9] = exp.(- ((delta_G_9) ./ kbT))
K[10] = exp.(- (delta_G_10 ./ kbT))
K[11] = exp.(- ((delta_G_11 ./ kbT)))
K[12] = exp.(- (delta_G_12 ./ kbT))
K[13] = exp.(- (delta_G_13 ./ kbT))

k_neg = k_pos ./ K

k_init = vcat(k_pos,k_neg)
end
