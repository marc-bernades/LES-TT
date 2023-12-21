function [rho_bw, mu_bw, T_tau_bw, T_bw, u_tau_bw,rho_tw, mu_tw, T_tau_tw, T_tw, u_tau_tw] = Calculate_wall_metrics(x,y,z,avg_rho,avg_mu,avg_T,avg_c_p,avg_kappa,avg_u,rho_b,T_b,u_b,num_points_x,num_points_y,num_points_z,delta)

% Make TKE symmetry on y
area_total_bw = 0;
rho_bw        = 0;
mu_bw         = 0;
T_bw          = 0;
c_p_bw        = 0;
kappa_bw      = 0;
u_boundary_bw = 0;
u_inner_bw    = 0;
T_boundary_bw = 0;
T_inner_bw    = 0;
area_total_tw = 0;
rho_tw        = 0;
mu_tw         = 0;
T_tw          = 0;
c_p_tw        = 0;
kappa_tw      = 0;
u_boundary_tw = 0;
u_inner_tw    = 0;
T_boundary_tw = 0;
T_inner_tw    = 0;


    for ii = 2:num_points_x-1
        for kk = 2:num_points_z-1

            % Bottom wall
            jj = 1;
            % Area bottom wall
            dx = 0.5*(x(ii+1,jj,kk) - x(ii-1,jj,kk));
            dz = 0.5*(z(ii,jj,kk+1) - z(ii,jj,kk-1));
            area_bw = dx.*dz;
            area_total_bw = area_total_bw + area_bw;
            % Parameters at bottom wall averaged inner and outer
            rho_bw   = rho_bw    + area_bw*0.5*(avg_rho(ii,jj,kk)   + avg_rho(ii,jj+1,kk));
            mu_bw    = mu_bw     + area_bw*0.5*(avg_mu(ii,jj,kk)    + avg_mu(ii,jj+1,kk));
            T_bw     = T_bw      + area_bw*0.5*(avg_T(ii,jj,kk)     + avg_T(ii,jj+1,kk));
            c_p_bw   = c_p_bw    + area_bw*0.5*(avg_c_p(ii,jj,kk)   + avg_c_p(ii,jj+1,kk));
            kappa_bw = kappa_bw  + area_bw*0.5*(avg_kappa(ii,jj,kk) + avg_kappa(ii,jj+1,kk));
            % Velocity and temperature inner and boundary
            u_boundary_bw  = u_boundary_bw + area_bw*avg_u(ii,jj,kk);
            u_inner_bw     = u_inner_bw    + area_bw*avg_u(ii,jj+1,kk);
            T_boundary_bw  = T_boundary_bw + area_bw*avg_T(ii,jj,kk);
            T_inner_bw     = T_inner_bw    + area_bw*avg_T(ii,jj+1,kk);


            % Top wall
            jj = num_points_y;
            % Area bottom wall
            dx = 0.5*(x(ii+1,jj,kk) - x(ii-1,jj,kk));
            dz = 0.5*(z(ii,jj,kk+1) - z(ii,jj,kk-1));
            area_tw = dx.*dz;
            area_total_tw = area_total_tw + area_tw;
            % Parameters at bottom wall averaged inner and outer
            rho_tw   = rho_tw    + area_tw*0.5*(avg_rho(ii,jj,kk)   + avg_rho(ii,jj-1,kk));
            mu_tw    = mu_tw     + area_tw*0.5*(avg_mu(ii,jj,kk)    + avg_mu(ii,jj-1,kk));
            T_tw     = T_tw      + area_tw*0.5*(avg_T(ii,jj,kk)     + avg_T(ii,jj-1,kk));
            c_p_tw   = c_p_tw    + area_tw*0.5*(avg_c_p(ii,jj,kk)   + avg_c_p(ii,jj-1,kk));
            kappa_tw = kappa_tw  + area_tw*0.5*(avg_kappa(ii,jj,kk) + avg_kappa(ii,jj-1,kk));
            % Velocity and temperature inner and boundary
            u_boundary_tw  = u_boundary_tw + area_tw*avg_u(ii,jj,kk);
            u_inner_tw     = u_inner_tw    + area_tw*avg_u(ii,jj-1,kk);
            T_boundary_tw  = T_boundary_tw + area_tw*avg_T(ii,jj,kk);
            T_inner_tw     = T_inner_tw    + area_tw*avg_T(ii,jj-1,kk);
        end

    end

    % Normalized by total area
    rho_bw        = rho_bw/area_total_bw;
    mu_bw         = mu_bw/area_total_bw;
    T_bw          = T_bw/area_total_bw;
    c_p_bw        = c_p_bw/area_total_bw;
    kappa_bw      = kappa_bw/area_total_bw;
    u_boundary_bw = u_boundary_bw/area_total_bw;
    u_inner_bw    = u_inner_bw/area_total_bw;
    T_boundary_bw = T_boundary_bw/area_total_bw;
    T_inner_bw    = T_inner_bw/area_total_bw;
    rho_tw        = rho_tw/area_total_tw;
    mu_tw         = mu_tw/area_total_tw;
    T_tw          = T_tw/area_total_tw;
    c_p_tw        = c_p_tw/area_total_tw;
    kappa_tw      = kappa_tw/area_total_tw;
    u_boundary_tw = u_boundary_tw/area_total_tw;
    u_inner_tw    = u_inner_tw/area_total_tw;
    T_boundary_tw = T_boundary_tw/area_total_tw;
    T_inner_tw    = T_inner_tw/area_total_tw;

    % Delta y
    delta_y_bw = y(1,2,1) - y(1,1,1);
    delta_y_tw = y(1,num_points_y,1) - y(1,num_points_y-1,1);

    % Tau wall
    tau_bw = mu_bw*( u_inner_bw - u_boundary_bw )/delta_y_bw;
    tau_tw = mu_tw*( u_inner_tw - u_boundary_tw )/delta_y_tw;

    %  u_tau
    u_tau_bw = sqrt( tau_bw/rho_bw );
    u_tau_tw = sqrt( tau_tw/rho_tw );

    % T_tau wall
    T_tau_bw = kappa_bw*( ( T_inner_bw - T_boundary_bw )/delta_y_bw )/( rho_bw*c_p_bw*u_tau_bw );
    T_tau_tw = kappa_tw*( ( T_boundary_tw - T_inner_tw )/delta_y_tw )/( rho_tw*c_p_tw*u_tau_tw );

    % alpha wall
    alpha_bw = kappa_bw/( rho_bw*c_p_bw );
    alpha_tw = kappa_tw/( rho_tw*c_p_tw );

    % Skin friction coefficient
    Cf_bw = tau_bw/( 0.5*rho_b*u_b*u_b );
    Cf_tw = tau_tw/( 0.5*rho_b*u_b*u_b );

    % Reynolds tau
    Re_tau_bw = rho_bw*u_tau_bw*delta/mu_bw;
    Re_tau_tw = rho_tw*u_tau_tw*delta/mu_tw;

    % Prandtl number at walls
    Pr_bw = c_p_bw*mu_bw/kappa_bw;
    Pr_tw = c_p_tw*mu_tw/kappa_tw;

    % Nusselt number at walls
    Nu_bw = delta*( ( T_boundary_bw - T_inner_bw )/delta_y_bw )/( T_bw - T_b );
    Nu_tw = delta*( ( T_boundary_tw - T_inner_tw )/delta_y_tw )/( T_tw - T_b );

    % Stanton number at walls
    St_bw = Nu_bw/( Re_tau_bw*Pr_bw );
    St_tw = Nu_tw/( Re_tau_tw*Pr_tw );




end