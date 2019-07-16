function kinetic_rate_mc_run(U,r)
try
    num_points=50


    U_vec = range(0.2,stop=1,length=num_points);
        #Calculation of Rate Constants
    k_pos = zeros(13)
    K = zeros(13)


        #Constant Constants
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
    G_U0s = Dict("O2aq"=>1.32+dGO2solv,"O2dl"=>1.32+dGO2solv,"O2"=>1.392+(1.2*dGOH),"OOH"=>1.25+dGOH,"O"=>-0.14+0.04+2*dGOH,"OH"=>-0.15+dGOH,"HOOH"=>1.533+0.04*dGOH,"Ob"=>0.5670+(2*dGOH),"OHb"=>0.1979+dGOH,"H2O2aq"=>1.75+H2O2aq_corr ," "=>0.,"b"=>0.)
    qs = Dict("O2aq"=>-4.,"O2dl"=>-4.,"O2"=>-4.,"OOH"=>-3.,"O"=>-2.,"OH"=>-1.,"HOOH"=>-2.,"Ob"=>-2.,"OHb"=>-1.,"H2O2aq"=>-2.," "=>0.,"b"=>0.)
    G_Us = Dict("O2aq"=>1.32+dGO2solv,"O2dl"=>1.32+dGO2solv,"O2"=>1.392+(1.2*dGOH),"OOH"=>1.25+dGOH,"O"=>-0.14+0.04+2*dGOH,"OH"=>-0.15+dGOH,"HOOH"=>1.533+0.04*dGOH,"Ob"=>0.5670+(2*dGOH),"OHb"=>0.1979+dGOH,"H2O2aq"=>1.75+H2O2aq_corr ," "=>0,"b"=>0.)
    U_vec = range(0.2,stop=1,length=num_points);
    rnd = zeros(12)
    rng=MersenneTwister(r)
    rnd = .05 .* randn(rng,12);
    #Potential: Input


    start_time = time_ns()
    U0=0.9

    deltaU = U-U0
    i=1
    for state in keys(G_U0s)
    	if(G_U0s[state]!=0)
    		G_Us[state] = G_U0s[state]+(qs[state]*deltaU)+rnd[i];
    	end
    	i+=1;
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
    #k_init = ones(26)
    x_o2aq = 2.34 .* (10 .^ -5)
    x_h2o = 1
    x_h2o2 = 0


        #Creation of parameter object
    p = p_base(k_init,x_o2aq,x_h2o,x_h2o2);

    #Problem setup & Initial Conditions
    diff_variables = trues(10);
    diff_variables[7]=false;
    diff_variables[10]=false

    tspan = (0,0.1)


    u_0 = [0.01,1.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]
    du_0 = [0.0,.1,.1,.1,.1,.1,0.1,.1,0.1,0.1]


    prob = DAEProblem(rate_equations,du_0,u_0,tspan,p,differential_vars=diff_variables)
    sol = solve(prob,IDA(max_num_iters_ic=10000))

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

    timetrack=(time_ns()-start_time)/1e9

        return o2dl,stara,o2stara,theta_oohstarA,theta_ostarA,theta_ohstarA,theta_h2o2starA,theta_ohstarB,theta_starB,theta_ostarB,timetrack

catch

    o2dl=-1.
    stara=-1.
    o2stara=-1.
    theta_oohstarA=-1.
    theta_ostarA=-1.
    theta_ohstarA=-1.
    theta_h2o2starA=-1.

    theta_ohstarB=-1.
    theta_starB=-1.
    theta_ostarB=-1.
    timetrack=-1.

    return o2dl,stara,o2stara,theta_oohstarA,theta_ostarA,theta_ohstarA,theta_h2o2starA,theta_ohstarB,theta_starB,theta_ostarB,timetrack
end




end
