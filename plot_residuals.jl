function plot_residuals(sol)
    time_vec = sol.t
    solution_vec = sol.u
    
    num_variables,num_timesteps=size(sol)
    resid_mat=zeros(num_timesteps-1,num_variables)
    
    for i = 2:num_timesteps
        for j = 1:num_variables
            resid_mat[i-1,j] = abs(sol[i][j]-sol[i-1][j])
            if(resid_mat[i-1,j]<1e-15)
                resid_mat[i-1,j]=1e-15
            end
        end
    end
    
    display(plot(sol.t[2:end],resid_mat,scale=:log10))
end
    
    