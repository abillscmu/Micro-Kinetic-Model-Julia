@everywhere using DifferentialEquations
@everywhere using Sundials
@everywhere using Random
@everywhere using DelimitedFiles
@everywhere using Distributions
@everywhere include("rate_equations.jl")
@everywhere include("datastructures.jl")
@everywhere include("kinetic_rate_mc_run.jl")
using SharedArrays

num_runs=1000
num_points=50
voltvec=range(.2,stop=1,length=num_points)
iterarr=[j for i=1:num_points, j=1:num_runs ]
voltarr=[voltvec[i] for i=1:num_points, j=1:num_runs]

o2dl,stara,o2stara,theta_oohstarA,theta_ostarA,theta_ohstarA,theta_h2o2starA,theta_ohstarB,theta_starB,theta_ostarB,timetrack=pmap(kinetic_rate_mc_run,voltarr,iterarr)



open("functiondir/o2dl.csv","w") do io
    writedlm(io,o2dl,',');
end

open("functiondir/stara.csv","w") do io
    writedlm(io,stara,',')
end

open("functiondir/o2stara.csv","w") do io
    writedlm(io,o2stara,',')
end

open("functiondir/theta_oohstarA.csv","w") do io
    writedlm(io,theta_oohstarA,',')
end

open("functiondir/theta_ostarA.csv","w") do io
    writedlm(io,theta_ostarA,',')
end

open("functiondir/theta_ohstarA.csv","w") do io
    writedlm(io,theta_ohstarA,',')
end

open("functiondir/theta_h2o2starA.csv","w") do io
    writedlm(io,theta_h2o2starA,',')
end

open("functiondir/theta_ohstarB.csv","w") do io
    writedlm(io,theta_ohstarB,',')
end

open("newdata/theta_starB.csv","w") do io
    writedlm(io,theta_starB,',')
end

open("functiondir/theta_ostarB.csv","w") do io
    writedlm(io,theta_ostarB,',')
end

open("functiondir/timetrack.csv","w") do io
   writedlm(io,timetrack,',')
end
