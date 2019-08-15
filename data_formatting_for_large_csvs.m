%% Data Formatting
%For some reason, matlab seems to not like large CSV's. 
num_voltage_points = 50;
num_monte_carlo = 100000;
[testrows,testcols]=size(o2dl);
if(testrows~=num_voltage_points||testcols~=num_monte_carlo)
    
    
    sizevec = [num_monte_carlo,num_voltage_points];
    o2dl = reshape(o2dl,sizevec)';
    o2stara=reshape(o2stara,sizevec)';
    stara=reshape(stara,sizevec)';
    theta_h2o2starA=reshape(theta_h2o2starA,sizevec)';
    theta_ohstarA=reshape(theta_ohstarA,sizevec)';
    theta_ohstarB=reshape(theta_ohstarB,sizevec)';
    theta_oohstarA=reshape(theta_oohstarA,sizevec)';
    theta_ostarA=reshape(theta_ostarA,sizevec)';
    theta_ostarB=reshape(theta_ostarB,sizevec)';
    theta_starB=reshape(theta_starB,sizevec)';
    t=reshape(t,sizevec)';
end
