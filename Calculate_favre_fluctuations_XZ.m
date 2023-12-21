function [R_favre_uu_fluct_calc,R_favre_vv_fluct_calc,R_favre_ww_fluct_calc] = Calculate_favre_fluctuations_XZ(u,v,w,avg_rho_xz,avg_rhou_xz,avg_rhov_xz,avg_rhow_xz,num_points_y)

R_favre_u_fluct_calc = zeros(size(u));
R_favre_v_fluct_calc = zeros(size(u));
R_favre_w_fluct_calc = zeros(size(u));

for jj = 2:num_points_y-1
    R_favre_u_fluct_calc(:,jj,:)  = u(:,jj,:) - avg_rhou_xz(jj)./avg_rho_xz(jj);
    R_favre_v_fluct_calc(:,jj,:)  = v(:,jj,:) - avg_rhov_xz(jj)./avg_rho_xz(jj);
    R_favre_w_fluct_calc(:,jj,:)  = w(:,jj,:) - avg_rhow_xz(jj)./avg_rho_xz(jj);
end
R_favre_uu_fluct_calc = R_favre_u_fluct_calc.^2;
R_favre_vv_fluct_calc = R_favre_v_fluct_calc.^2;
R_favre_ww_fluct_calc = R_favre_w_fluct_calc.^2;

end