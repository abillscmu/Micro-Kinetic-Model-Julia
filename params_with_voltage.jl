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

    ##NEED U HERE (THIS WILL BE THE BEGINING OF THE LOOP!)
U=0.2

    deltaU = U-U0
    for state in keys(G_U0s)
        G_Us[state] = G_U0s[state]+(qs[state]*deltaU);
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
    s2 = exp.(-beta * 0.2244)
    s3 = 1.
    s4 = 1.
    s13 = s2


    G_0_3 = -0.9+delta_G_3
    G_0_4 = -0.9+delta_G_4
    G_0_5 = -0.9+delta_G_5
    G_0_6 = -0.9+delta_G_6
    G_0_8 = -0.9+delta_G_8
    G_0_9 = -0.9+delta_G_9
    G_0_11 = -0.9+delta_G_11

    E_1 = 0
    E_2 = 0
    E_7 = 0.48
    E_10 = 0.37
    E_12 = 0.46
    E_13 = 0

    #Delta G error 0.05
    #lognormal exponential of Gaussian

    U_3 = -G_0_3
    U_4 = -G_0_4
    U_5 = -G_0_5
    U_6 = -G_0_6
    U_8 = -G_0_8
    U_9 = -G_0_9
    U_11 = -G_0_11
    A = 1.0 .* (10 .^ 9)



    #Potential: Input


        k_pos[1] = (8 .* (10 .^ 5))
        k_pos[2] = kbTh *f0 *s2 * min(1.,exp(-beta*delta_G_2))
        k_pos[3] = kbTh * f0 * s3 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_3 ) ))
        k_pos[4] = kbTh * f0 * s4 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_4 ) ))
        k_pos[5] = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_5 ) ))
        k_pos[6] = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_6 ) ))
        k_pos[7] = kbTh * min(1.,exp( - beta * ( 0.48 + 0.69 * ( 0.9255 + delta_G_7) ) ))
        k_pos[8] = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_8 ) ))
        k_pos[9] = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_9 ) )))
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
        #for n = 1:26
        #        if abs(k_init[n]>1e20)
        #            k_init[n]=1e20
        #        end
        #    end

        #Creation of parameter object
        p = p_base(k_init,x_o2aq,x_h2o,x_h2o2);
