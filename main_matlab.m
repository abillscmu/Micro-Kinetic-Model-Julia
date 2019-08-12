%%Main Kinetic Rate Function
function out = main_matlab()

num_points=50;
o2dl=zeros(num_points,1);
stara=zeros(num_points,1);
o2stara=zeros(num_points,1);
theta_oohstarA=zeros(num_points,1);
theta_ostarA=zeros(num_points,1);
theta_ohstarA=zeros(num_points,1);
theta_h2o2starA=zeros(num_points,1);
theta_ohstarB=zeros(num_points,1);
theta_starB=zeros(num_points,1);
theta_ostarB=zeros(num_points,1);



U_vec = linspace(0.2,1,num_points);
h=4.14e-15;
%Calculation of Rate Constants
k_pos = zeros(13,1);
K = zeros(13,1);

%Constant Constants
kb = 8.617 .* (10 .^ -5);
T = 298;
kbT = kb .* T;
kbTh = kbT ./ h;

%Trib
E = 0.26;
beta = 0.5;

cH2O = 1e3./(16+((2*1.008)));
kH_H2O=1.3e-3;
dGO2solv=kbT*log(cH2O/kH_H2O);
dGH2O2solv=kbT*log(cH2O);
dGOH=0;
U0=0.9;
OHB_destabilization=0.36;
Ob_destabilization=0.0;
H2O2aq_corr=dGH2O2solv;
O2ads_corr=0;

elements={'O2aq','O2dl','O2','OOH','O','OH','HOOH','Ob','OHb','H2O2aq',' ','b'};
G_U0s_values = {1.32+dGO2solv,1.32+dGO2solv,1.392+(1.2*dGOH),1.25+dGOH,-0.14+0.04+2*dGOH,-0.15+dGOH,1.533+0.04*dGOH,0.5670+(2*dGOH),0.1979+dGOH,1.75+H2O2aq_corr,0,0,};
G_U0s=containers.Map(elements,G_U0s_values);
qs_values={-4,-4,-4,-3,-2,-1,-2,-2,-1,-2,0,0};
qs = containers.Map(elements,qs_values);
G_Us = containers.Map(elements,G_U0s_values);



%Potential: Input
for n = 1:num_points
    U = U_vec(n);
    deltaU = U-U0;
    for i=1:length(G_Us)
        G_Us(elements{i}) = G_U0s(elements{i})+(qs(elements{i})*deltaU);
    end
    
    delta_G_1 = G_Us('O2dl')-G_Us('O2aq');
    delta_G_2 = G_Us('O2')-G_Us('O2dl');
    delta_G_3 = G_Us('OOH')-G_Us('O2');
    delta_G_4 = G_Us('O')-G_Us('OOH');
    delta_G_5 = G_Us('OH')-G_Us('O');
    delta_G_6 = G_Us(' ')-G_Us('OH');
    delta_G_7 = G_Us('O')+G_Us('Ob')-G_Us('O2');
    delta_G_8 = G_Us('OHb')-G_Us('Ob');
    delta_G_9 = G_Us('b')-G_Us('OHb');
    delta_G_10 = G_Us('Ob')+G_Us('OH')-G_Us('OOH');
    delta_G_11 = G_Us('HOOH')-G_Us('OOH');
    delta_G_12 = G_Us('OH')+G_Us('OHb')-G_Us('HOOH');
    delta_G_13 = G_Us('H2O2aq')-G_Us('HOOH');
    
    beta = 1. / ( kbT );
    
    f0 = exp(-beta * 0.2244);
    f1 = 1;
    s2 = exp(-beta * 0.05916);
    s3 = 1.;
    s4 = 1.;
    s13 = s2;
    
    k_pos(1) = (8 .* (10 .^ 5));
    k_pos(2) = kbTh *f0 *s2 * min(1.,exp(-beta*delta_G_2));
    k_pos(3) = kbTh * f0 * s3 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_3 ) ));
    k_pos(4) = kbTh * f0 * s4 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_4 ) ));
    k_pos(5) = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_5 ) ));
    k_pos(6) = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_6 ) ));
    k_pos(7) = kbTh * min(1.,exp( - beta * ( 0.48 + 0.69 * ( 0.9255 + delta_G_7) ) ));
    k_pos(8) = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_8 ) ));
    k_pos(9) = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_9 ) ));
    k_pos(10) = kbTh * min(1.,exp( - beta * ( 0.37 + 0.39 * (0.8330 + delta_G_10) ) ));
    k_pos(11) = kbTh * f0 * min(1.,exp( - beta * ( 0.26 + 0.5 * delta_G_11 ) ));
    k_pos(12) =  kbTh * min(1.,exp( - beta * ( 0.462 + 0.19 * (1.4853 + delta_G_12) ) ));
    k_pos(13) = kbTh * s13 * f0 * min(1.,exp( - beta * ( delta_G_13 ) ));
    
    K(1) = exp(- (delta_G_1 ./ kbT));
    K(2) = exp(- (delta_G_2 ./ kbT));
    K(3) = exp(- ((delta_G_3) ./ kbT));
    K(4) = exp(- ((delta_G_4) ./ kbT));
    K(5) = exp(- ((delta_G_5) ./ kbT));
    K(6) = exp(- ((delta_G_6) ./ kbT));
    K(7) = exp(- (delta_G_7 ./ kbT));
    K(8) = exp(- ((delta_G_8) ./ kbT));
    K(9) = exp(- ((delta_G_9) ./ kbT));
    K(10) = exp(- (delta_G_10 ./ kbT));
    K(11) = exp(- ((delta_G_11 ./ kbT)));
    K(12) = exp(- (delta_G_12 ./ kbT));
    K(13) = exp(- (delta_G_13 ./ kbT));
    
    k_neg = k_pos ./ K;
    k_init = vertcat(k_pos,k_neg);
    x_o2aq = 2.34 .* (10 .^ -5);
    x_h2o = 1;
    x_h2o2 = 0;
    
    %Creation of parameter object
    p=[k_init;x_o2aq;x_h2o;x_h2o2];
    
    %Problem setup & Initial Conditions
    diff_variables = ones(10,1);
    diff_variables(7)=0;
    diff_variables(10)=0;
    M = diag(diff_variables);
    
    %MATLAB DAE Interface
    options = odeset('Mass',M);
    tspan = [0 10];
    %println(k_init(11))
    %println(k_init(24))
    u_0 = [0.01,1.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]';
    prob = @(t,y)rate_equations(t,y,p);
    tic
    [t,y] = ode15s(prob,tspan,u_0,options);
    toc
    
    o2dl(n) = y(end,1);
    stara(n) = y(end,2);
    o2stara(n) = y(end,7);
    theta_oohstarA(n) = y(end,4);
    theta_ostarA(n) = y(end,5);
    theta_ohstarA(n) = y(end,6);
    theta_h2o2starA(n) = y(end,3);
    theta_ohstarB(n) = y(end,8);
    theta_starB(n) = y(end,10);
    theta_ostarB(n) = y(end,9);
    
    
    %Solve and extract outputs
    
end
%{
figure(1)
clf
semilogy(U_vec,o2dl)
hold on
semilogy(U_vec,stara);
semilogy(U_vec,o2stara);
semilogy(U_vec,theta_oohstarA);
semilogy(U_vec,theta_ostarA)
semilogy(U_vec,theta_ohstarA);
semilogy(U_vec,theta_h2o2starA);
legend({'o2dl','stara','o2stara','theta_oohstara','theta_ostarA','theta_ohstara','theta_h2o2stara'})


%}
out=0;
end

