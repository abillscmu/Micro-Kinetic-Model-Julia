using Random
function mc_modifier(prob,i,something)
    h=prob.p.h
    #Calculation of Rate Constants
    k_pos = zeros(13)
    K = zeros(13)

    #  Constant Constants
    kb = prob.p.kb
    T = prob.p.T
    kbT = kb .* T
    kbTh = kbT ./ h

    #Trib
    E = prob.p.E
    beta = prob.p.beta

    #Create Randomness
    #rng = MersenneTwister(Int(round(time()*i)))
    rnd = zeros(13)


    delta_G_1 =  prob.p.delta_G[1] + rnd[1]
    delta_G_2 =  prob.p.delta_G[2] + rnd[2]
    delta_G_3 =  prob.p.delta_G[3] + rnd[3]
    delta_G_4 =  prob.p.delta_G[4] + rnd[4]
    delta_G_5 =  prob.p.delta_G[5] + rnd[5]
    delta_G_6 =  prob.p.delta_G[6] + rnd[6]
    delta_G_7 =  prob.p.delta_G[7] + rnd[7]
    delta_G_8 =  prob.p.delta_G[8] + rnd[8]
    delta_G_9 =  prob.p.delta_G[9] + rnd[9]
    delta_G_10 =  prob.p.delta_G[10] + rnd[10]
    delta_G_11 =  prob.p.delta_G[11] + rnd[11]
    delta_G_12 =  prob.p.delta_G[12] + rnd[12]
    delta_G_13 =  prob.p.delta_G[13] + rnd[13]


    G_0_3 = -0.9+delta_G_3
    G_0_4 = -0.9+delta_G_4
    G_0_5 = -0.9+delta_G_5
    G_0_6 = -0.9+delta_G_6
    G_0_8 = -0.9+delta_G_8
    G_0_9 = -0.9+delta_G_9
    G_0_11 = -0.9+delta_G_11

    E_1 = prob.p.chem_E[1]
    E_2 = prob.p.chem_E[2]
    E_7 = prob.p.chem_E[3]
    E_10 = prob.p.chem_E[4]
    E_12 = prob.p.chem_E[5]
    E_13 = prob.p.chem_E[6]

    #Delta G error 0.05
    #lognormal exponential of Gaussian

    U_3 = -G_0_3
    U_4 = -G_0_4
    U_5 = -G_0_5
    U_6 = -G_0_6
    U_8 = -G_0_8
    U_9 = -G_0_9
    U_11 = -G_0_11
    A = prob.p.A


    U = prob.p.U;

    k_pos[1] = (8 .* (10 .^ 5)) .* exp.(-E_1 ./ kbT )
    k_pos[2] = (1 .* (10 .^ 8)) .* exp.(-E_2 ./ kbT )
    k_pos[3] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_3))) ./ kbT)
    k_pos[4] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_4))) ./ kbT)
    k_pos[5] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_5))) ./ kbT)
    k_pos[6] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_6))) ./ kbT)
    k_pos[7] = kbTh .* exp.(-E_7 ./ kbT)
    k_pos[8] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_8))) ./ kbT)
    k_pos[9] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_9))) ./ kbT)
    k_pos[10] = kbTh .* exp.(-E_10 ./ kbT)
    k_pos[11] = A .* exp.(- (E ./ kbT)) .* exp.(- ((beta .* (U-U_11))) ./ kbT)
    k_pos[12] = kbTh .* exp.(-E_12./ kbT)
    k_pos[13] = (1.00 .* (10 .^ 8)) .* exp.(-E_13 ./ kbT)

    K[1] = exp.(- (delta_G_1 ./ kbT))
    K[2] = exp.(- (delta_G_2 ./ kbT))
    K[3] = exp.(- ((G_0_3+U) ./ kbT))
    K[4] = exp.(- ((G_0_4+U) ./ kbT))
    K[5] = exp.(- ((G_0_5+U) ./ kbT))
    K[6] = exp.(- ((G_0_6+U) ./ kbT))
    K[7] = exp.(- (delta_G_3 ./ kbT))
    K[8] = exp.(- ((G_0_8+U) ./ kbT))
    K[9] = exp.(- ((G_0_9+U) ./ kbT))
    K[10] = exp.(- (delta_G_10 ./ kbT))
    K[11] = exp.(- ((G_0_11+U) ./ kbT))
    K[12] = exp.(- (delta_G_12 ./ kbT))
    K[13] = exp.(- (delta_G_13 ./ kbT))

    k_neg = k_pos ./ K

    k_init = vcat(k_pos,k_neg)
    prob.p.k=k_init
    return prob

end
