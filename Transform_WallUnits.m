function [y_plus_bw,u_plus_bw,T_plus_bw,y_plus_tw,u_plus_tw,T_plus_tw] = Transform_WallUnits(y,avg_u,avg_T,u_tau_bw,rho_bw,mu_bw,T_bw,T_tau_bw,u_tau_tw,rho_tw,mu_tw,T_tw,T_tau_tw,num_points_x,num_points_y,num_points_z,delta)

num_points_xz     = num_points_x*num_points_z;

% Make TKE symmetry on y
y_plus_bw        = zeros(1,num_points_y);
u_plus_bw        = zeros(1,num_points_y);
T_plus_bw        = zeros(1,num_points_y);
y_plus_tw        = zeros(1,num_points_y);
u_plus_tw        = zeros(1,num_points_y);
T_plus_tw        = zeros(1,num_points_y);


for jj = 2:num_points_y-1
    for ii = 2:num_points_x-1
        for kk = 2:num_points_z-1

            % Top wall
            aux_j = num_points_y - jj;
            y_plus_tw(aux_j) = y_plus_tw(aux_j) + (2*delta - y(ii,jj,kk))*(u_tau_tw/(mu_tw/rho_tw))/num_points_xz;
            u_plus_tw(aux_j) = u_plus_tw(aux_j) + avg_u(ii,jj,kk)/u_tau_tw/num_points_xz;
            T_plus_tw(aux_j) = T_plus_tw(aux_j) + (avg_T(ii,jj,kk) - T_tw)/T_tau_tw/num_points_xz;

            % Bottom wall
            aux_j = jj;
            y_plus_bw(aux_j) = y_plus_bw(aux_j) + y(ii,jj,kk)*(u_tau_bw/(mu_bw/rho_bw))/num_points_xz;
            u_plus_bw(aux_j) = u_plus_bw(aux_j) + avg_u(ii,jj,kk)/u_tau_bw/num_points_xz;
            T_plus_bw(aux_j) = T_plus_bw(aux_j) + (avg_T(ii,jj,kk) - T_bw)/T_tau_bw/num_points_xz;

        
        end


    end
end




end