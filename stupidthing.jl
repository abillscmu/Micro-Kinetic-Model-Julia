
indices = Dict("O2aq"=>1,"O2dl"=>2,"O2"=>3,"OOH"=>4,"O"=>5,"OH"=>6,"HOOH"=>7,"Ob"=>8,"OHb"=>9,"H2O2aq"=>10 ," "=>11,"b"=>12)
values2  = Dict("O2aq"=>1.32+dGO2solv,"O2dl"=>1.32+dGO2solv,"O2"=>1.392+(1.2*dGOH),"OOH"=>1.25+dGOH,"O"=>-0.14+0.04+2*dGOH,"OH"=>-0.15+dGOH,"HOOH"=>1.533+0.04*dGOH,"Ob"=>0.5670+(2*dGOH),"OHb"=>0.1979+dGOH,"H2O2aq"=>1.75+H2O2aq_corr ," "=>0.,"b"=>0.)
qs = Dict("O2aq"=>-4.,"O2dl"=>-4.,"O2"=>-4.,"OOH"=>-3.,"O"=>-2.,"OH"=>-1.,"HOOH"=>-2.,"Ob"=>-2.,"OHb"=>-1.,"H2O2aq"=>-2.," "=>0.,"b"=>0.)
###ORDERING


for state in keys(indices)
    println(G_U0s[state]==params3[indices[state]+1])

end
