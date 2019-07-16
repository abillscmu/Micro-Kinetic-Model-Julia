function plottinginterval(matrix,color,fignum,legendname)
[num_voltage_points,num_monte_carlo_runs]=size(matrix);






for n = 1:num_voltage_points
    vec = matrix(n,:);
    vec=vec(vec~=0);
    svec = sort(vec);
    li = round(0.05.*length(vec));
    ui=round(length(vec))-li;
    low(n)=svec(li);
    high(n)=svec(ui);
    ave(n) = mean(svec(li:ui));
    
    
        
end
legendentry = [legendname ' average'];
legendci = [legendname ' 95\% CI'];
xvec=linspace(0.2,1,num_voltage_points);
figure(fignum)
ci = [color '--'];
semilogy(xvec,low,ci,'LineWidth',2)
hold on
semilogy(xvec,high,ci,'LineWidth',2)
semilogy(xvec,ave,color,'LineWidth',3,'DisplayName',legendentry)
end