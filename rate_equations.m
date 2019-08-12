%Rate Equations
%Alec Bills, June 10
%Source: Hansen et. al 2014
function out = rate_equations(t,u,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%DECOMPOSE CONSTANTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %k's
    k_1=k(1);
    k_2=k(2);
    k_3=k(3);
    k_4=k(4);
    k_5=k(5);
    k_6=k(6);
    k_7=k(7);
    k_8=k(8);
    k_9=k(9);
    k_10=k(10);
    k_11=k(11);
    k_12=k(12);
    k_13=k(13);
    k_m1=k(14);
    k_m2=k(15);
    k_m3=k(16);
    k_m4=k(17);
    k_m5=k(18);
    k_m6=k(19);
    k_m7=k(20);
    k_m8=k(21);
    k_m9=k(22);
    k_m10=k(23);
    k_m11=k(24);
    k_m12=k(25);
    k_m13=k(26);
    
    %x's
    x_o2aq=k(27);
    x_h2o=k(28);
    x_h2o2=k(29);
    
%%%%%%%%%%%%%%%%%%%%%%%%%DECOMPOSE PARAMETERS IN THE SOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_o2dl=u(1);
    theta_starA=u(2);
    theta_o2starA=u(7);
    theta_oohstarA=u(4);
    theta_ostarA=u(5);
    theta_ohstarA=u(6);
    theta_h2o2starA=u(3);
    
    theta_ohstarB=u(8);
    theta_starB=u(10);
    theta_ostarB=u(9);
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%CONSTRUCT EQUATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dx_o2dldt = ((k_1) .* (x_o2aq)) - ((k_m1) .* (x_o2dl)) - ((k_2) .* (x_o2dl) .* (theta_starA)) + ((k_m2 .* theta_o2starA)); %Equation 1

    dtheta_starAdt = ((k_6) .* (theta_ohstarA)) - ((k_m6) .* (theta_starA) .* (x_h2o)) - ((k_2) .* (x_o2dl) .* (theta_starA)) + ((k_m2) .* (theta_o2starA)) + ((k_13) .* (theta_h2o2starA)) - ((k_m13) .* (theta_starA) .* (x_h2o2)); %Equation 2

    dtheta_o2starAdt = ((k_2) .* (x_o2dl) .* (theta_starA)) - ((k_m2) .* (theta_o2starA)) - ((k_3) .* (theta_o2starA)) + ((k_m3) .* (theta_oohstarA)) - ((k_7) .* (theta_o2starA) .* (theta_starB)) + ((k_m7) .* (theta_ostarA) .* (theta_ostarB));%Equation 3 %REPLACED WITH CONSERVATION EQUATION

    dtheta_oohstarAdt = ((k_3).*(theta_o2starA)) - ((k_m3) .* (theta_oohstarA)) - ((k_4) .* (theta_oohstarA)) + ((k_m4) .* (x_h2o) .* (theta_ostarA)) - ((k_10) .* (theta_oohstarA) .* (theta_starB)) + ((k_m10) .* (theta_ohstarA) .* (theta_ostarB)) - ((k_11) .* (theta_oohstarA)) + ((k_m11) .* (theta_h2o2starA));%Equation 4

    dtheta_ostarAdt = ((k_4) .* (theta_oohstarA)) - ((k_m4) .* (x_h2o) .* (theta_ostarA)) - ((k_5) .* (theta_ostarA)) + ((k_m5) .* (theta_ohstarA)) + ((k_7) .* (theta_o2starA) .* (theta_starB)) - ((k_m7) .* (theta_ostarA) .* (theta_ostarB));%Equation 5

    dtheta_ohstarAdt = ((k_5) .* (theta_ostarA)) - ((k_m5) .* (theta_ohstarA)) - ((k_6) .* (theta_ohstarA)) + ((k_m6) .* (theta_starA) .* (x_h2o)) + ((k_10) .* (theta_oohstarA) .* (theta_starB)) - ((k_m10) .* (theta_ohstarA) .* (theta_ostarB)) + ((k_12) .* (theta_h2o2starA) .* (theta_starB)) - ((k_m12) .* (theta_ohstarA) .* (theta_ohstarB)); %Equation 6

    dtheta_h2o2starAdt = ((k_11) .* (theta_oohstarA)) - ((k_m11) .* (theta_h2o2starA)) - ((k_12) .* (theta_h2o2starA) .* (theta_starB)) + ((k_m12) .* (theta_ohstarA) .* (theta_ohstarB)) - ((k_13) .* (theta_h2o2starA)) + ((k_m13) .* (theta_starA) .* (x_h2o2));%Equation 7

    dtheta_ostarBdt = ((k_7) .* (theta_o2starA) .* (theta_starB)) - ((k_m7) .* (theta_ostarA) .* (theta_ostarB)) + ((k_10) .* (theta_oohstarA) .* (theta_starB)) - ((k_m10) .* (theta_ohstarA) .* (theta_ostarB)) - ((k_8) .* (theta_ostarB)) + ((k_m8) .* (theta_ohstarB));%Equation 8

    dtheta_ohstarBdt = ((k_8) .* (theta_ostarB)) - ((k_m8) .* (theta_ohstarB)) + ((k_12) .* (theta_h2o2starA) .* (theta_starB)) - ((k_m12) .* (theta_ohstarA) .* (theta_ohstarB)) - ((k_9) .* (theta_ohstarB)) + ((k_m9) .* (x_h2o) .* (theta_starB)); %Equation 9

    dtheta_starBdt = ((k_9) .* (theta_ohstarB)) - ((k_m9) .* (x_h2o) .* (theta_starB)) - ((k_7) .*(theta_o2starA) .* (theta_starB)) + ((k_m7) .* (theta_ostarA) .* (theta_ostarB)) - ((k_10) .* (theta_oohstarA) .* (theta_starB)) + ((k_m10) .* (theta_ohstarA) .* (theta_ostarB)) - ((k_12) .* (theta_h2o2starA) .* (theta_starB)) + ((k_m12) .* (theta_ohstarA) .* (theta_ohstarB));%Equation 10
    
    
%%%%%%%%%%%%%%%%%%%%%%%OUTPUT EQUATIONS IN USABLE FORM (CONSTRUCT DIFFERENTIAL OUTPUT VECTOR) %%%%%%%%%%%%%%%%%%%%%%%%
    out = zeros(10,1);
    out(1)=dx_o2dldt;
    out(2)=dtheta_starAdt;
    %out(3)=dtheta_o2starAdt-du(3)
    out(7)=u(2)+u(3)+u(4)+u(5)+u(6)+u(7)-1.0;
    out(4)=dtheta_oohstarAdt;
    out(5)=dtheta_ostarAdt;
    out(6)=dtheta_ohstarAdt;
    out(3)=dtheta_h2o2starAdt;
    out(8)=dtheta_ohstarBdt;
    out(9)=dtheta_ostarBdt;
    out(10)=u(8)+u(9)+u(10)-1.0;

end