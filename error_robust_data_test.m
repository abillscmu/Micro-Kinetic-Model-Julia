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






%% Find where things are zero
o2dlmask = o2dl==0;
o2staramask = o2stara==0;
staramask=stara==0;
h2o2staramask=theta_h2o2starA==0;
ohstaramask=theta_ohstarA==0;
ohstarBmask=theta_ohstarB==0;
oohstaramask=theta_oohstarA==0;
ostaramask=theta_ostarA==0;
ostarbmask=theta_ostarB==0;
starbmask=theta_starB==0;
tmask = t==0;

allmat=o2dlmask+o2staramask+staramask+h2o2staramask+ohstaramask+ohstarBmask+oohstaramask+ostaramask+ostarbmask+starbmask+tmask;
allmask=allmat==11;
anymask=allmat>0;

%%
non_broken_t = t(~anymask);
outliers_removed_t=non_broken_t(non_broken_t<1);



