function [avg_rho,avg_mu,avg_u,avg_v,avg_w,avg_rhou,avg_rhov,avg_rhow] = Spatial_avg_XZ(rho,mu,u,v,w,num_points_x,num_points_y,num_points_z)

num_points_xz     = (num_points_x-2)*(num_points_z-2);

% Make TKE symmetry on y
avg_rho        = zeros(1,num_points_y);
avg_mu         = zeros(1,num_points_y);
avg_u          = zeros(1,num_points_y);
avg_v          = zeros(1,num_points_y);
avg_w          = zeros(1,num_points_y);
avg_rhou       = zeros(1,num_points_y);
avg_rhov       = zeros(1,num_points_y);
avg_rhow       = zeros(1,num_points_y);



for jj = 2:num_points_y-1
    for ii = 2:num_points_x-1
        for kk = 2:num_points_z-1

            % Compute across entire channel
            aux_j = jj;
            avg_rho(aux_j)        = avg_rho(aux_j)  + rho(ii,jj,kk)/num_points_xz;
            avg_mu(aux_j)        = avg_mu(aux_j)    + mu(ii,jj,kk)/num_points_xz;
            avg_u(aux_j)          = avg_u(aux_j)    + u(ii,jj,kk)/num_points_xz;
            avg_v(aux_j)          = avg_v(aux_j)    + v(ii,jj,kk)/num_points_xz;
            avg_w(aux_j)          = avg_w(aux_j)    + w(ii,jj,kk)/num_points_xz;
            avg_rhou(aux_j)       = avg_rhou(aux_j) + rho(ii,jj,kk).*u(ii,jj,kk)/num_points_xz;
            avg_rhov(aux_j)       = avg_rhov(aux_j) + rho(ii,jj,kk).*v(ii,jj,kk)/num_points_xz;
            avg_rhow(aux_j)       = avg_rhow(aux_j) + rho(ii,jj,kk).*w(ii,jj,kk)/num_points_xz;
           
        end

    end
end


end