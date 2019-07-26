%% Monte Carlo UQ Data Visualization
filepath = 'newdata/';

%% Read Data
o2dl=csvread([filepath 'o2dl.csv']);
o2stara=csvread([filepath 'o2stara.csv']);
stara=csvread([filepath 'stara.csv']);
theta_h2o2starA=csvread([filepath 'theta_h2o2starA.csv']);
theta_ohstarA=csvread([filepath 'theta_ohstarA.csv']);
theta_ohstarB=csvread([filepath 'theta_ohstarB.csv']);
theta_oohstarA=csvread([filepath 'theta_oohstarA.csv']);
theta_ostarA=csvread([filepath 'theta_ostarA.csv']);
theta_ostarB=csvread([filepath 'theta_ostarB.csv']);
theta_starB=csvread([filepath 'theta_starB.csv']);
t=csvread([filepath 'timetrack.csv']);

%% Plot Data
figure(1)
clf
plottinginterval(o2dl,'g',1,'o2dl')
plottinginterval(o2stara,'b',1,'o2stara')
plottinginterval(stara,'r',1,'stara')
plottinginterval(theta_oohstarA,'c',1,'oohstara')
plottinginterval(theta_ostarA,'m',1,'ostara')
plottinginterval(theta_ohstarA,'y',1,'ohstara')
plottinginterval(theta_h2o2starA,'k',1,'h2o2stara')
title('Monte Carlo UQ');
ylabel('Theta')
xlabel('U')
set(gca,'FontSize',20)