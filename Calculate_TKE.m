function [TKE, avg_TKE_y, avg_TKE_tw, avg_TKE_bw, avg_R_favre_uu_bw, avg_R_favre_vv_bw, avg_R_favre_ww_bw, ...
    avg_R_favre_uu_tw, avg_R_favre_vv_tw, avg_R_favre_ww_tw] = Calculate_TKE(R_favre_uu,R_favre_vv,R_favre_ww,num_points_x,num_points_y,num_points_z)

num_points_xz     = num_points_x*num_points_z;

% Make TKE symmetry on y

avg_R_favre_uu_tw = zeros(1,ceil(num_points_y));
avg_R_favre_vv_tw = zeros(1,ceil(num_points_y));
avg_R_favre_ww_tw = zeros(1,ceil(num_points_y));
avg_R_favre_uu_bw = zeros(1,ceil(num_points_y));
avg_R_favre_vv_bw = zeros(1,ceil(num_points_y));
avg_R_favre_ww_bw = zeros(1,ceil(num_points_y));
avg_R_favre_uu_y  = zeros(1,num_points_y);
avg_R_favre_vv_y  = zeros(1,num_points_y);
avg_R_favre_ww_y  = zeros(1,num_points_y);

avg_TKE_y   = zeros(1,num_points_y);
avg_TKE_tw  = zeros(1,num_points_y);
avg_TKE_bw  = zeros(1,num_points_y);


% TKE point by point
TKE    = 1/2*(R_favre_uu + R_favre_vv + R_favre_ww);

% Idx to where calculate TKE
idx = 1;

for jj = (idx+1):(num_points_y-idx)
    for ii = (idx+1):(num_points_x-idx)
        for kk = (idx+1):(num_points_z-idx)

            % Compute across entire channel
            aux_j = jj;
            avg_R_favre_uu_y(aux_j) = avg_R_favre_uu_y(aux_j) + R_favre_uu(ii,jj,kk)/num_points_xz;
            avg_R_favre_vv_y(aux_j) = avg_R_favre_vv_y(aux_j) + R_favre_vv(ii,jj,kk)/num_points_xz;
            avg_R_favre_ww_y(aux_j) = avg_R_favre_ww_y(aux_j) + R_favre_ww(ii,jj,kk)/num_points_xz;
            avg_TKE_y(aux_j)        = avg_TKE_y(aux_j)        + TKE(ii,jj,kk)/num_points_xz;


            % Top wall
            %             if jj > ceil(0.5*num_points_y)
            aux_j = num_points_y - jj + 1;
            avg_R_favre_uu_tw(aux_j) = avg_R_favre_uu_tw(aux_j) + R_favre_uu(ii,jj,kk)/num_points_xz;
            avg_R_favre_vv_tw(aux_j) = avg_R_favre_vv_tw(aux_j) + R_favre_vv(ii,jj,kk)/num_points_xz;
            avg_R_favre_ww_tw(aux_j) = avg_R_favre_ww_tw(aux_j) + R_favre_ww(ii,jj,kk)/num_points_xz;
            avg_TKE_tw(aux_j)        = avg_TKE_tw(aux_j)        + TKE(ii,jj,kk)/num_points_xz;

            %             else
            % Bottom wall
            aux_j  = jj;
            avg_R_favre_uu_bw(aux_j) = avg_R_favre_uu_bw(aux_j) + R_favre_uu(ii,jj,kk)/num_points_xz;
            avg_R_favre_vv_bw(aux_j) = avg_R_favre_vv_bw(aux_j) + R_favre_vv(ii,jj,kk)/num_points_xz;
            avg_R_favre_ww_bw(aux_j) = avg_R_favre_ww_bw(aux_j) + R_favre_ww(ii,jj,kk)/num_points_xz;
            avg_TKE_bw(aux_j)        = avg_TKE_bw(aux_j)        + TKE(ii,jj,kk)/num_points_xz;

            %             end

        end

    end
end


% TKE across all y
TKE_y  = 1/2*(avg_R_favre_uu_y + avg_R_favre_vv_y + avg_R_favre_ww_y);

% TKE for each wall
TKE_tw = 1/2*(avg_R_favre_uu_tw + avg_R_favre_vv_tw + avg_R_favre_ww_tw);
TKE_bw = 1/2*(avg_R_favre_uu_bw + avg_R_favre_vv_bw + avg_R_favre_ww_bw);


end