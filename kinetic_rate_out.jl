
function calc_ans(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13)
input_vec=vcat(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13);

k_init=calc_ks(input_vec)
#k_init = ones(26)
x_o2aq = 2.34 .* (10 .^ -5)
x_h2o = 1
x_h2o2 = 0


#Creation of parameter object
#Parameter Structure: U,g
p = vcat(k_init,x_o2aq,x_h2o,x_h2o2);



tspan = (0,10.)

u_0 = [1.,0.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]
du_0 = [0.0,.1,-.1,.1,.1,.1,0.1,.1,0.1,0.1]
mm = ones(10);
mm[7]=0
mm[10]=0
MM=Diagonal(mm);
func = ODEFunction(rate_equations,mass_matrix=MM)
prob = ODEProblem(func,u_0,tspan,p)
sol = solve(prob,Rodas4())

o2dl=sol.u[end][1]
stara=sol.u[end][2]
o2stara=sol.u[end][7]
theta_oohstarA=sol.u[end][4]
theta_ostarA=sol.u[end][5]
theta_ohstarA=sol.u[end][6]
theta_h2o2starA=sol.u[end][3]

theta_ohstarB=sol.u[end][8]
theta_starB=sol.u[end][10]
theta_ostarB=sol.u[end][9]
return [o2dl,stara,o2stara,theta_oohstarA,theta_ostarA,theta_ohstarA,theta_h2o2starA,theta_ohstarB,theta_starB,theta_ostarB]
end

function calc_ans(input_vec::Array{Float32,1})

k_init=calc_ks(input_vec)
#k_init = ones(26)
x_o2aq = 2.34 .* (10 .^ -5)
x_h2o = 1
x_h2o2 = 0


#Creation of parameter object
#Parameter Structure: U,g
p = vcat(k_init,x_o2aq,x_h2o,x_h2o2);



tspan = (0,10.)

u_0 = [1.,0.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]
du_0 = [0.0,.1,-.1,.1,.1,.1,0.1,.1,0.1,0.1]
mm = ones(10);
mm[7]=0
mm[10]=0
MM=Diagonal(mm);
func = ODEFunction(rate_equations,mass_matrix=MM)
prob = ODEProblem(func,u_0,tspan,p)
sol = solve(prob,Rodas4())

o2dl=sol.u[end][1]
stara=sol.u[end][2]
o2stara=sol.u[end][7]
theta_oohstarA=sol.u[end][4]
theta_ostarA=sol.u[end][5]
theta_ohstarA=sol.u[end][6]
theta_h2o2starA=sol.u[end][3]

theta_ohstarB=sol.u[end][8]
theta_starB=sol.u[end][10]
theta_ostarB=sol.u[end][9]
return [o2dl,stara,o2stara,theta_oohstarA,theta_ostarA,theta_ohstarA,theta_h2o2starA,theta_ohstarB,theta_starB,theta_ostarB]
end

function calc_ans(input_vec::Array{Float64,1})

k_init=calc_ks(input_vec)
#k_init = ones(26)
x_o2aq = 2.34 .* (10 .^ -5)
x_h2o = 1
x_h2o2 = 0


#Creation of parameter object
#Parameter Structure: U,g
p = vcat(k_init,x_o2aq,x_h2o,x_h2o2);



tspan = (0,10.)

u_0 = [1.,0.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]
du_0 = [0.0,.1,-.1,.1,.1,.1,0.1,.1,0.1,0.1]
mm = ones(10);
mm[7]=0
mm[10]=0
MM=Diagonal(mm);
func = ODEFunction(rate_equations,mass_matrix=MM)
prob = ODEProblem(func,u_0,tspan,p)
sol = solve(prob,Rodas4())

o2dl=sol.u[end][1]
stara=sol.u[end][2]
o2stara=sol.u[end][7]
theta_oohstarA=sol.u[end][4]
theta_ostarA=sol.u[end][5]
theta_ohstarA=sol.u[end][6]
theta_h2o2starA=sol.u[end][3]

theta_ohstarB=sol.u[end][8]
theta_starB=sol.u[end][10]
theta_ostarB=sol.u[end][9]
return [o2dl,stara,o2stara,theta_oohstarA,theta_ostarA,theta_ohstarA,theta_h2o2starA,theta_ohstarB,theta_starB,theta_ostarB]
end
