%%Main Kinetic Rate Function 
function main_matlab()

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
E = 0.26 + 0.05;
beta = 0.5;

delta_G_1 = 0;
delta_G_2 = -0.2 + 0.05;
delta_G_3 = -.14 + 0.05;
delta_G_4 = -1.35 + 0.05;
delta_G_5 = -0.05 + 0.05;
delta_G_6 = 0.15 + 0.05;
delta_G_7 = -.93 + 0.05;
delta_G_8 = -.38 + 0.05;
delta_G_9 = -0.2 + 0.05;
delta_G_10 = -0.83 + 0.05;
delta_G_11 = 0.28 + 0.05;
delta_G_12 = -1.49 + 0.05;
delta_G_13 = 0.32 + 0.05;


G_0_3 = -0.9+delta_G_3;
G_0_4 = -0.9+delta_G_4;
G_0_5 = -0.9+delta_G_5;
G_0_6 = -0.9+delta_G_6;
G_0_8 = -0.9+delta_G_8;
G_0_9 = -0.9+delta_G_9;
G_0_11 = -0.9+delta_G_11;

E_1 = 0;
E_2 = 0;
E_7 = 0.48 + 0.05;
E_10 = 0.37 + 0.05;
E_12 = 0.46 + 0.05;
E_13 = 0;



U_3 = -G_0_3;
U_4 = -G_0_4;
U_5 = -G_0_5;
U_6 = -G_0_6;
U_8 = -G_0_8;
U_9 = -G_0_9;
U_11 = -G_0_11;
A = 1.0 .* (10 .^ 9);



%Potential: Input
for n = 1:num_points
    U = U_vec(n);    
    k_pos(1) = (8 .* (10 .^ 5)) .* exp(-E_1 ./ kbT );
    k_pos(2) = (1 .* (10 .^ 8)) .* exp(-E_2 ./ kbT );    
    k_pos(3) = A .* exp(- (E ./ kbT)) .* exp(- ((beta .* (U-U_3))) ./ kbT);
    k_pos(4) = A .* exp(- (E ./ kbT)) .* exp(- ((beta .* (U-U_4))) ./ kbT);    
    k_pos(5) = A .* exp(- (E ./ kbT)) .* exp(- ((beta .* (U-U_5))) ./ kbT);
    k_pos(6) = A .* exp(- (E ./ kbT)) .* exp(- ((beta .* (U-U_6))) ./ kbT);   
    k_pos(7) = kbTh .* exp(-E_7 ./ kbT);
    k_pos(8) = A .* exp(- (E ./ kbT)) .* exp(- ((beta .* (U-U_8))) ./ kbT);   
    k_pos(9) = A .* exp(- (E ./ kbT)) .* exp(- ((beta .* (U-U_9))) ./ kbT);
    k_pos(10) = kbTh .* exp(-E_10 ./ kbT);    
    k_pos(11) = A .* exp(- (E ./ kbT)) .* exp(- ((beta .* (U-U_11))) ./ kbT);
    k_pos(12) = kbTh .* exp(-E_12./ kbT);    
    k_pos(13) = (1.00 .* (10 .^ 8)) .* exp(-E_13 ./ kbT);
    K(1) = exp(- (delta_G_1 ./ kbT));
    K(2) = exp(- (delta_G_2 ./ kbT));    
    K(3) = exp(- ((G_0_3+U) ./ kbT));
    K(4) = exp(- ((G_0_4+U) ./ kbT));    
    K(5) = exp(- ((G_0_5+U) ./ kbT));
    K(6) = exp(- ((G_0_6+U) ./ kbT));    
    K(7) = exp(- (delta_G_3 ./ kbT));
    K(8) = exp(- ((G_0_8+U) ./ kbT));    
    K(9) = exp(- ((G_0_9+U) ./ kbT));
    K(10) = exp(- (delta_G_10 ./ kbT));    
    K(11) = exp(- ((G_0_11+U) ./ kbT));
    K(12) = exp(- (delta_G_12 ./ kbT));    
    K(13) = exp(- (delta_G_13 ./ kbT));
    
    k_neg = k_pos ./ K;    
    k_init = vertcat(k_pos,k_neg);
    x_o2aq = 2.34 .* (10 .^ -5);    
    x_h2o = 1;
    x_h2o2 = 0;
    
    %Creation of parameter object
    p.k = k_init;    
    p.x_o2aq = x_o2aq;
    p.x_h2o = x_h2o;    
    p.x_h2o2 = x_h2o2;
    
    %Problem setup & Initial Conditions
    diff_variables = ones(10,1);    
    diff_variables(7)=0;
    diff_variables(10)=0;   
    M = diag(diff_variables);
    
    %MATLAB DAE Interface
    options = odeset('Mass',M);
    tspan = [0 0.1];
    %println(k_init(11))
    %println(k_init(24))
    u_0 = [0.01,1.,0.0,0.0,0.0,0.0,0.0,1.,0.0,0.0]';
    prob = @(t,y)rate_equations(t,y,p);
    [t,y] = ode15s(prob,tspan,u_0,options);
    
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
semilogy(U_vec,o2dl);
hold on
semilogy(U_vec,stara);
semilogy(U_vec,o2stara);
semilogy(U_vec,theta_oohstarA);
semilogy(U_vec,theta_ostarA);
semilogy(U_vec,theta_ohstarA);
semilogy(U_vec,theta_h2o2starA);
legend({'o2dl','stara','o2stara','oohstara','ostarA','ohstarA','h2o2starA'});
%}

