%% Model implementation
% Fi defined by unclosed_terms source run code
delta_filt  = delta*fi; % Fi defined in unclosed_terms as governs the _filt fields
Filter_type = 'CDLF_Box';
HP_model    = 'HighPressure';

% Delta x, y and z based on CentralDerivative_d1_2ndOrder
[~,dx,~] = CentralDerivative_d1_2ndOrder(x);
[dy,~,~] = CentralDerivative_d1_2ndOrder(y);
[~,~,dz] = CentralDerivative_d1_2ndOrder(z);

% Aditional filter terms to compute C_I
rho_filt_FG          = FilterFields(rho_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhouu_filt_favre_FG  = FilterFields(rho_filt.*u_filt_favre.*u_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhouv_filt_favre_FG  = FilterFields(rho_filt.*u_filt_favre.*v_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhouw_filt_favre_FG  = FilterFields(rho_filt.*u_filt_favre.*w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhovv_filt_favre_FG  = FilterFields(rho_filt.*v_filt_favre.*v_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhovw_filt_favre_FG  = FilterFields(rho_filt.*v_filt_favre.*w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhoww_filt_favre_FG  = FilterFields(rho_filt.*w_filt_favre.*w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);

rhou_filt_favre_FG   = FilterFields(rho_filt.*u_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhov_filt_favre_FG   = FilterFields(rho_filt.*v_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhow_filt_favre_FG   = FilterFields(rho_filt.*w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);

[du_y_filt_favre,du_x_filt_favre,du_z_filt_favre] = CentralDerivative_d1_2ndOrder(u_filt_favre);
[dv_y_filt_favre,dv_x_filt_favre,dv_z_filt_favre] = CentralDerivative_d1_2ndOrder(v_filt_favre);
[dw_y_filt_favre,dw_x_filt_favre,dw_z_filt_favre] = CentralDerivative_d1_2ndOrder(w_filt_favre);

% Strain tensor
divU_filt_favre = du_x_filt_favre./dx + dv_y_filt_favre./dy + dw_z_filt_favre./dz;
S_uu_filt_favre = (du_x_filt_favre./dx + du_x_filt_favre./dx) - 2/3*divU_filt_favre;
S_uv_filt_favre = (du_y_filt_favre./dy + dv_x_filt_favre./dx);
S_uw_filt_favre = (du_z_filt_favre./dz + dw_x_filt_favre./dx);
S_vv_filt_favre = (dv_y_filt_favre./dy + dv_y_filt_favre./dy) - 2/3*divU_filt_favre;
S_vw_filt_favre = (dv_z_filt_favre./dz + dw_y_filt_favre./dy);
S_ww_filt_favre = (dw_z_filt_favre./dz + dw_z_filt_favre./dz) - 2/3*divU_filt_favre;

% Absolute strain tensor
S_ij_filt_favre_abs = (1/2*(S_uu_filt_favre.^2 + 2*S_uv_filt_favre.^2 + 2*S_uw_filt_favre.^2 + S_vv_filt_favre.^2 + 2*S_vw_filt_favre.^2 + S_ww_filt_favre.^2)).^0.5;

% FG level filtered
S_ij_filt_favre_abs_FG            = FilterFields(S_ij_filt_favre_abs,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
S_ij_filt_favre_abs_2_rho_filt_FG = FilterFields(rho_filt.*S_ij_filt_favre_abs.^2,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhouu_kk_filt_favre_FG            = FilterFields(rho_filt.*(u_filt_favre.^2 + v_filt_favre.^2 + w_filt_favre.^2),2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhou_kk_filt_favre_FG             = FilterFields(rho_filt.*(u_filt_favre + v_filt_favre + w_filt_favre),2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);


%% Trace model Yoshizawa
% C_I_num = (rhouu_kk_filt_favre_FG - 1./rho_filt_FG.*rhou_kk_filt_favre_FG.^2);
% C_I_den = (rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_abs_FG.^2 - delta_filt.^2.*S_ij_filt_favre_abs_2_rho_filt_FG);

C_I_x_num = (rhouu_filt_favre_FG - 1./rho_filt_FG.*rhou_filt_favre_FG.^2);
C_I_x_den = (rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_abs_FG.^2 - delta_filt.^2.*S_ij_filt_favre_abs_2_rho_filt_FG);
C_I_y_num = (rhovv_filt_favre_FG - 1./rho_filt_FG.*rhov_filt_favre_FG.^2);
C_I_y_den = (rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_abs_FG.^2 - delta_filt.^2.*S_ij_filt_favre_abs_2_rho_filt_FG);
C_I_z_num = (rhoww_filt_favre_FG - 1./rho_filt_FG.*rhow_filt_favre_FG.^2);
C_I_z_den = (rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_abs_FG.^2 - delta_filt.^2.*S_ij_filt_favre_abs_2_rho_filt_FG);

% [C_I_num_avg] = Spatial_avg_XZ_var(C_I_num,num_points_x,num_points_y,num_points_z, fi, idx_filt);
% [C_I_den_avg] = Spatial_avg_XZ_var(C_I_den,num_points_x,num_points_y,num_points_z, fi, idx_filt);

[C_I_x_num_avg] = Spatial_avg_XZ_var(C_I_x_num,num_points_x,num_points_y,num_points_z, fi, idx_filt);
[C_I_x_den_avg] = Spatial_avg_XZ_var(C_I_x_den,num_points_x,num_points_y,num_points_z, fi, idx_filt);
[C_I_y_num_avg] = Spatial_avg_XZ_var(C_I_y_num,num_points_x,num_points_y,num_points_z, fi, idx_filt);
[C_I_y_den_avg] = Spatial_avg_XZ_var(C_I_y_den,num_points_x,num_points_y,num_points_z, fi, idx_filt);
[C_I_z_num_avg] = Spatial_avg_XZ_var(C_I_z_num,num_points_x,num_points_y,num_points_z, fi, idx_filt);
[C_I_z_den_avg] = Spatial_avg_XZ_var(C_I_z_den,num_points_x,num_points_y,num_points_z, fi, idx_filt);

% C_I_avg = C_I_num_avg./C_I_den_avg; C_I_avg(isnan(C_I_avg)) = 0;

C_I_x_avg = C_I_x_num_avg./C_I_x_den_avg;
C_I_y_avg = C_I_y_num_avg./C_I_y_den_avg;
C_I_z_avg = C_I_z_num_avg./C_I_z_den_avg;


figure;
subplot(3,1,1)
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, C_I_x_avg(idx_filt:end-(idx_filt)));
xlabel('$y/\delta$','interpreter','latex')
ylabel('${C_I}_x$','interpreter','latex')
subplot(3,1,2)
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, C_I_y_avg(idx_filt:end-(idx_filt)));
xlim([0 2])
xlabel('$y/\delta$','interpreter','latex')
ylabel('${C_I}_y$','interpreter','latex')
subplot(3,1,3)
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, C_I_z_avg(idx_filt:end-(idx_filt)));
xlabel('$y/\delta$','interpreter','latex')
ylabel('${C_I}_z$','interpreter','latex')
xlim([0 2])

% Initialize trace
% Tau_kk   = zeros(size(x));
Tau_kk_x = zeros(size(x));
Tau_kk_y = zeros(size(x));
Tau_kk_z = zeros(size(x));

for jj = 1:length(Tau_kk_x(2,:,2))
%     Tau_kk(:,jj,:)   = C_I_avg(jj).*delta_filt(:,jj,:).^2.*S_ij_filt_favre_abs(:,jj,:).^2;
    Tau_kk_x(:,jj,:) = C_I_x_avg(jj).*delta_filt(:,jj,:).^2.*S_ij_filt_favre_abs(:,jj,:).^2;
    Tau_kk_y(:,jj,:) = C_I_y_avg(jj).*delta_filt(:,jj,:).^2.*S_ij_filt_favre_abs(:,jj,:).^2;
    Tau_kk_z(:,jj,:) = C_I_z_avg(jj).*delta_filt(:,jj,:).^2.*S_ij_filt_favre_abs(:,jj,:).^2;
end

Tau_kk_u_SFS{1} = Tau_kk_x;
Tau_kk_v_SFS{1} = Tau_kk_y;
Tau_kk_w_SFS{1} = Tau_kk_z;
Tau_kk_SFS{1}   = Tau_kk_u_SFS{1} + Tau_kk_v_SFS{1} + Tau_kk_w_SFS{1};

% Tau_kk_u_SFS{1} = Tau_kk/3;
% Tau_kk_v_SFS{1} = Tau_kk/3;
% Tau_kk_w_SFS{1} = Tau_kk/3;
% Tau_kk_SFS{1}   = Tau_kk;


%% Vreman model
C_w          = 0.325; % Homogeneous isotropic turbulence Nicoud and Ducros 1999 for WALE model
S_u2_sum     = S_uu_filt_favre.^2 + S_uv_filt_favre.^2 + S_uw_filt_favre.^2;
S_v2_sum     = S_uv_filt_favre.^2 + S_vv_filt_favre.^2 + S_vw_filt_favre.^2;
S_w2_sum     = S_uw_filt_favre.^2 + S_vw_filt_favre.^2 + S_ww_filt_favre.^2;

Tau_kk_u_SFS{2} = C_w*delta_filt.^2.*S_u2_sum;
Tau_kk_v_SFS{2} = C_w*delta_filt.^2.*S_v2_sum;
Tau_kk_w_SFS{2} = C_w*delta_filt.^2.*S_w2_sum;
Tau_kk_SFS{2}   = Tau_kk_u_SFS{2} + Tau_kk_v_SFS{2} + Tau_kk_w_SFS{2};

% Tau_kk_SFS{2}   = C_w*delta_filt.^2.*(S_u2_sum + S_v2_sum + S_w2_sum);
% Tau_kk_u_SFS{2} = Tau_kk_SFS{2}/3;
% Tau_kk_v_SFS{2} = Tau_kk_SFS{2}/3;
% Tau_kk_w_SFS{2} = Tau_kk_SFS{2}/3;

%% Load DNS filt Tau (Turbulent Stress Tensor)
Tau_kk_DNS_filt = zeros(size(x));
Tau_xx_DNS_filt = Tau_xz; Tau_xy_DNS_filt = Tau_xy; Tau_xz_DNS_filt = Tau_xz;
Tau_yy_DNS_filt = Tau_yy; Tau_yz_DNS_filt = Tau_yz; Tau_zz_DNS_filt = Tau_zz;

for ii = idx_filt:num_points_x-(idx_filt-1)
    for jj = idx_filt:num_points_y-(idx_filt-1) % Only filtered points from y_outer
        for kk = idx_filt:num_points_z-(idx_filt-1)

            tau_filt_ij = [Tau_xx(ii,jj,kk), Tau_xy(ii,jj,kk), Tau_xz(ii,jj,kk);
                Tau_xy(ii,jj,kk), Tau_yy(ii,jj,kk), Tau_yz(ii,jj,kk);
                Tau_xz(ii,jj,kk), Tau_yz(ii,jj,kk), Tau_zz(ii,jj,kk)];

            tau_filt_kk = tau_filt_ij(1,1) + tau_filt_ij(2,2) + tau_filt_ij(3,3);

            Tau_xx_DNS_filt(ii,jj,kk) = Tau_xx(ii,jj,kk) - 1/3*tau_filt_kk;
            Tau_yy_DNS_filt(ii,jj,kk) = Tau_yy(ii,jj,kk) - 1/3*tau_filt_kk;
            Tau_zz_DNS_filt(ii,jj,kk) = Tau_zz(ii,jj,kk) - 1/3*tau_filt_kk;

            Tau_kk_DNS_filt(ii,jj,kk) = tau_filt_kk;

        end
    end
end
            

figure; hold on
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_kk_DNS_filt(idx_filt+30,idx_filt:end-(idx_filt),idx_filt+30))
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_kk_SFS{1}(idx_filt+30,idx_filt:end-(idx_filt),idx_filt+30))
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_kk_SFS{2}(idx_filt+30,idx_filt:end-(idx_filt),idx_filt+30))
legend('DNS filt', 'Yoshizawa','Vreman')
xlabel('$y/\delta$','interpreter','latex')
ylabel('${\tau}_{kk}$','interpreter','latex')
xlim([0 2])

%% Smagorinsky model
C_s = 0.1; % Turbulent channel flow Deadorff 1970
ve_Smagorinsky = C_s^2.*delta_filt.^2.*S_ij_filt_favre_abs;


Tau_xx_SFS{1} = -ve_Smagorinsky.*S_uu_filt_favre;
Tau_xy_SFS{1} = -ve_Smagorinsky.*S_uv_filt_favre;
Tau_xz_SFS{1} = -ve_Smagorinsky.*S_uw_filt_favre;
Tau_yy_SFS{1} = -ve_Smagorinsky.*S_vv_filt_favre;
Tau_yz_SFS{1} = -ve_Smagorinsky.*S_vw_filt_favre;
Tau_zz_SFS{1} = -ve_Smagorinsky.*S_ww_filt_favre;


%% Dynamic model
u_filt_favre_FG   = FilterFields(rho_filt.*u_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_FG;
v_filt_favre_FG   = FilterFields(rho_filt.*v_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_FG;
w_filt_favre_FG   = FilterFields(rho_filt.*w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_FG;

% Compute strain tensor at FG level
[du_y_filt_favre_FG,du_x_filt_favre_FG,du_z_filt_favre_FG] = CentralDerivative_d1_2ndOrder(u_filt_favre_FG);
[dv_y_filt_favre_FG,dv_x_filt_favre_FG,dv_z_filt_favre_FG] = CentralDerivative_d1_2ndOrder(v_filt_favre_FG);
[dw_y_filt_favre_FG,dw_x_filt_favre_FG,dw_z_filt_favre_FG] = CentralDerivative_d1_2ndOrder(w_filt_favre_FG);

divU_filt_favre_FG = du_x_filt_favre_FG./dx + dv_y_filt_favre_FG./dy + dw_z_filt_favre_FG./dz;
S_uu_filt_favre_FG = (du_x_filt_favre_FG./dx + du_x_filt_favre_FG./dx) - 2/3*divU_filt_favre_FG;
S_uv_filt_favre_FG = (du_y_filt_favre_FG./dy + dv_x_filt_favre_FG./dx);
S_uw_filt_favre_FG = (du_z_filt_favre_FG./dz + dw_x_filt_favre_FG./dx);
S_vv_filt_favre_FG = (dv_y_filt_favre_FG./dy + dv_y_filt_favre_FG./dy) - 2/3*divU_filt_favre_FG;
S_vw_filt_favre_FG = (dv_z_filt_favre_FG./dz + dw_y_filt_favre_FG./dy);
S_ww_filt_favre_FG = (dw_z_filt_favre_FG./dz + dw_z_filt_favre_FG./dz) - 2/3*divU_filt_favre_FG;

% Absolute strain tensor at FG level
S_ij_filt_favre_FG_abs = (1/2*(S_uu_filt_favre_FG.^2 + 2*S_uv_filt_favre_FG.^2 + 2*S_uw_filt_favre_FG.^2 + S_vv_filt_favre_FG.^2 + 2*S_vw_filt_favre_FG.^2 + S_ww_filt_favre_FG.^2)).^0.5;


M_uu_FG = FilterFields(rho_filt.*delta_filt.^2.*S_ij_filt_favre_abs.*S_uu_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
M_uu    = -rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_FG_abs.*S_uu_filt_favre_FG + M_uu_FG;
M_uv_FG = FilterFields(rho_filt.*delta_filt.^2.*S_ij_filt_favre_abs.*S_uv_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
M_uv    = -rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_FG_abs.*S_uv_filt_favre_FG + M_uv_FG;
M_uw_FG = FilterFields(rho_filt.*delta_filt.^2.*S_ij_filt_favre_abs.*S_uw_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
M_uw    = -rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_FG_abs.*S_uw_filt_favre_FG + M_uw_FG;
M_vv_FG = FilterFields(rho_filt.*delta_filt.^2.*S_ij_filt_favre_abs.*S_vv_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
M_vv    = -rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_FG_abs.*S_vv_filt_favre_FG + M_vv_FG;
M_vw_FG = FilterFields(rho_filt.*delta_filt.^2.*S_ij_filt_favre_abs.*S_vw_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
M_vw    = -rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_FG_abs.*S_vw_filt_favre_FG + M_vw_FG;
M_ww_FG = FilterFields(rho_filt.*delta_filt.^2.*S_ij_filt_favre_abs.*S_ww_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
M_ww    = -rho_filt_FG.*(2*delta_filt).^2.*S_ij_filt_favre_FG_abs.*S_ww_filt_favre_FG + M_ww_FG;

rhouu_filt_FG     = FilterFields(rho_filt.*u_filt_favre.*u_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhouv_filt_FG     = FilterFields(rho_filt.*u_filt_favre.*v_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhouw_filt_FG     = FilterFields(rho_filt.*u_filt_favre.*w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhovv_filt_FG     = FilterFields(rho_filt.*v_filt_favre.*v_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhovw_filt_FG     = FilterFields(rho_filt.*v_filt_favre.*w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhoww_filt_FG     = FilterFields(rho_filt.*w_filt_favre.*w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);

rhou_filt        = rho_filt.*u_filt_favre;
rhov_filt        = rho_filt.*v_filt_favre;
rhow_filt        = rho_filt.*w_filt_favre;

rhou_filt_FG   = FilterFields(rhou_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhov_filt_FG   = FilterFields(rhov_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rhow_filt_FG   = FilterFields(rhow_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);


L_uu = rhouu_filt_favre_FG - rhou_filt_FG.*rhou_filt_FG./rho_filt_FG;
L_uv = rhouv_filt_favre_FG - rhou_filt_FG.*rhov_filt_FG./rho_filt_FG;
L_uw = rhouw_filt_favre_FG - rhou_filt_FG.*rhow_filt_FG./rho_filt_FG;
L_vv = rhovv_filt_favre_FG - rhov_filt_FG.*rhov_filt_FG./rho_filt_FG;
L_vw = rhovw_filt_favre_FG - rhov_filt_FG.*rhow_filt_FG./rho_filt_FG;
L_ww = rhoww_filt_favre_FG - rhow_filt_FG.*rhow_filt_FG./rho_filt_FG;

% % C_d component wise (tensor 3x3) and least squares
% C_d_uu = L_uu.*M_uu./(M_uu.*M_uu);
% C_d_uv = L_uv.*M_uv./(M_uv.*M_uv);
% C_d_uw = L_uw.*M_uw./(M_uw.*M_uw);
% C_d_vv = L_vv.*M_vv./(M_vv.*M_vv);
% C_d_vw = L_vw.*M_vw./(M_vw.*M_vw);
% C_d_ww = L_ww.*M_ww./(M_ww.*M_ww);

% % Plot evolution dynamic coefficient
% figure; hold on
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, C_d_uu(idx_filt:end-(idx_filt)))
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, C_d_uv(idx_filt:end-(idx_filt)))
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, C_d_uw(idx_filt:end-(idx_filt)))
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, C_d_vv(idx_filt:end-(idx_filt)))
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, C_d_vw(idx_filt:end-(idx_filt)))
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, C_d_ww(idx_filt:end-(idx_filt)))
% xlim([0 2])
% legend('uu','uv','uw','vv','vw','ww')
% xlabel('$y/\delta$','interpreter','latex')
% ylabel('${C_d}$','interpreter','latex')

% % C_d minimum square error satisfying each grid point
% C_d = zeros(size(rho_filt));
% 
% for ii = 1:length(C_d_uu(:,1,1))
%     for jj = 1:length(C_d_uu(1,:,1))
%         for kk = 1:length(C_d_uu(1,1,:))
%             C_d(ii,jj,kk) = mean([C_d_uu(ii,jj,kk),C_d_vv(ii,jj,kk),C_d_ww(ii,jj,kk),C_d_uv(ii,jj,kk),C_d_uv(ii,jj,kk),C_d_uw(ii,jj,kk),C_d_uw(ii,jj,kk),C_d_vw(ii,jj,kk),C_d_vw(ii,jj,kk)]);
%         end
%     end
% end

% C_d inner product
C_d = (M_uu.*L_uu + 2*M_uv.*L_uv + 2*M_uw.*L_uw + M_vv.*L_vv + 2*M_vw.*L_vw + M_ww.*L_ww)./(M_uu.*M_uu + 2*M_uv.*M_uv + 2*M_uw.*M_uw + M_vv.*M_vv + 2*M_vw.*M_vw + M_ww.*M_ww);
C_d(C_d<0) = 0;


% Ensemble average
C_d_avg = Spatial_avg_XZ_var(C_d,num_points_x,num_points_y,num_points_z, fi, idx_filt);
C_d_avg(C_d_avg < 0) = 0;


% Dynamic viscosity
ve_Smagorinsky_dyn = zeros(size(x));
for jj = 1:length(C_d_avg)
    ve_Smagorinsky_dyn(:,jj,:) = C_d_avg(jj)^2.*delta_filt(:,jj,:).^2.*S_ij_filt_favre_abs(:,jj,:);
end

% Turbulent stress tensor
Tau_xx_SFS{2} = -ve_Smagorinsky_dyn.*S_uu_filt_favre;
Tau_xy_SFS{2} = -ve_Smagorinsky_dyn.*S_uv_filt_favre;
Tau_xz_SFS{2} = -ve_Smagorinsky_dyn.*S_uw_filt_favre;
Tau_yy_SFS{2} = -ve_Smagorinsky_dyn.*S_vv_filt_favre;
Tau_yz_SFS{2} = -ve_Smagorinsky_dyn.*S_vw_filt_favre;
Tau_zz_SFS{2} = -ve_Smagorinsky_dyn.*S_ww_filt_favre;

%% Anisotropic minimum-dissipation model (AMD)
Omega_uu_filt_favre = (du_x_filt_favre./dx - du_x_filt_favre./dx);
Omega_uv_filt_favre = (du_y_filt_favre./dy - dv_x_filt_favre./dx);
Omega_uw_filt_favre = (du_z_filt_favre./dz - dw_x_filt_favre./dx);
Omega_vu_filt_favre = -Omega_uv_filt_favre;
Omega_vv_filt_favre = (dv_y_filt_favre./dy - dv_y_filt_favre./dy);
Omega_vw_filt_favre = (dv_z_filt_favre./dz - dw_y_filt_favre./dy);
Omega_wu_filt_favre = -Omega_uw_filt_favre;
Omega_wv_filt_favre = -Omega_uv_filt_favre;
Omega_ww_filt_favre = (dw_z_filt_favre./dz - dw_z_filt_favre./dz);

ve_AMR = zeros(size(x));

% Decomposition point by point
for ii = idx_filt:num_points_x-(idx_filt-1)
    for jj = idx_filt:num_points_y-(idx_filt-1) % Only filtered points from y_outer
        for kk = idx_filt:num_points_z-(idx_filt-1)

            S_ij = [S_uu_filt_favre(ii,jj,kk), S_uv_filt_favre(ii,jj,kk), S_uw_filt_favre(ii,jj,kk);
                S_uv_filt_favre(ii,jj,kk), S_vv_filt_favre(ii,jj,kk), S_vw_filt_favre(ii,jj,kk);
                S_uw_filt_favre(ii,jj,kk), S_vw_filt_favre(ii,jj,kk), S_ww_filt_favre(ii,jj,kk)];

            Omega_ij = [Omega_uu_filt_favre(ii,jj,kk), Omega_uv_filt_favre(ii,jj,kk), Omega_uw_filt_favre(ii,jj,kk);
                Omega_vu_filt_favre(ii,jj,kk), Omega_vv_filt_favre(ii,jj,kk), Omega_vw_filt_favre(ii,jj,kk);
                Omega_wu_filt_favre(ii,jj,kk), Omega_wv_filt_favre(ii,jj,kk), Omega_ww_filt_favre(ii,jj,kk)];

            S_ij_2           = 1/2*S_ij*1/2*S_ij;
            S_ij_3           = 1/2*S_ij*1/2*S_ij*1/2*S_ij;
            Omega_ij_2       = 1/2*Omega_ij*1/2*Omega_ij;
            S_ij_Omega_ij_2  = 1/2*S_ij*Omega_ij_2;

            % Trace
            S_kk_2 = S_ij_2(1,1) + S_ij_2(2,2) + S_ij_2(3,3);
            S_kk_3 = S_ij_3(1,1) + S_ij_3(2,2) + S_ij_3(3,3);

            Omega_kk_2   = Omega_ij_2(1,1) + Omega_ij_2(2,2) + Omega_ij_2(3,3);
            S_Omega_kk_2 = S_ij_Omega_ij_2(1,1) + S_ij_Omega_ij_2(2,2) + S_ij_Omega_ij_2(3,3);


            % Tensor invariants
            I_1 = S_kk_2;
            I_2 = Omega_kk_2;
            I_3 = S_kk_3;
            I_4 = S_Omega_kk_2;

            % Viscosity model
            ve_AMR(ii,jj,kk) = max(0,-(I_3 - I_4))./(I_1 - I_2);
        end
    end
end



% C_delta = 1./(4./(fi.*dx.^2) + 4./(fi.*dy.^2) + 4./(fi.*dz.^2));
C_delta = 0.3; % Impose C = 0.3 to give sufficient dissipation and use filte width
ve_AMR  = ve_AMR.*(C_delta.*delta_filt).^2;


Tau_xx_SFS{3} = - ve_AMR.*S_uu_filt_favre;
Tau_xy_SFS{3} = - ve_AMR.*S_uv_filt_favre;
Tau_xz_SFS{3} = - ve_AMR.*S_uw_filt_favre;
Tau_yy_SFS{3} = - ve_AMR.*S_vv_filt_favre;
Tau_yz_SFS{3} = - ve_AMR.*S_vw_filt_favre;
Tau_zz_SFS{3} = - ve_AMR.*S_ww_filt_favre;


%% WALE Model
% Velocity gradient
A_uu_filt_favre = du_x_filt_favre./dx; A_uv_filt_favre = du_y_filt_favre./dy; A_uw_filt_favre = du_z_filt_favre./dz;
A_vu_filt_favre = dv_x_filt_favre./dx; A_vv_filt_favre = dv_y_filt_favre./dy; A_vw_filt_favre = dv_z_filt_favre./dz;
A_wu_filt_favre = dw_x_filt_favre./dx; A_wv_filt_favre = dw_y_filt_favre./dy; A_ww_filt_favre = dw_z_filt_favre./dz;


% Grid point based
% Initialize D_Sigma operator
D_WALE         = zeros(size(dx));
C_s_wale_2_num = zeros(size(dx));
C_s_wale_2_den = zeros(size(dx));

for ii = idx_filt:num_points_x-(idx_filt-1)
    for jj = idx_filt:num_points_y-(idx_filt-1) % Only filtered points from y_outer
        for kk = idx_filt:num_points_z-(idx_filt-1)

            S_ij = 0.5*[S_uu_filt_favre(ii,jj,kk) S_uv_filt_favre(ii,jj,kk) S_uw_filt_favre(ii,jj,kk);
                S_uv_filt_favre(ii,jj,kk) S_vv_filt_favre(ii,jj,kk) S_vw_filt_favre(ii,jj,kk);
                S_uw_filt_favre(ii,jj,kk) S_vw_filt_favre(ii,jj,kk) S_ww_filt_favre(ii,jj,kk)];

            A_ij = [A_uu_filt_favre(ii,jj,kk) A_uv_filt_favre(ii,jj,kk) A_uw_filt_favre(ii,jj,kk);
                A_uv_filt_favre(ii,jj,kk) A_vv_filt_favre(ii,jj,kk) A_vw_filt_favre(ii,jj,kk);
                A_uw_filt_favre(ii,jj,kk) A_vw_filt_favre(ii,jj,kk) A_ww_filt_favre(ii,jj,kk)];

            A_ji = A_ij';

            S_ij_d = 0.5*((A_ij*A_ij + A_ji*A_ji) - 2/3*trace(A_ij*A_ij)*eye(3));

            % Inner product S_ij_d S_ij_d
            S_ij_d_2 = sum(sum(S_ij_d.*S_ij_d)); 

            % Inner product S_ij S_ij
            S_ij_2  = sum(sum(S_ij.*S_ij)); 

            % WALE Operator
            D_WALE(ii,jj,kk) = (S_ij_d_2.^(3/2))./(S_ij_2.^(5/2) + S_ij_d_2.^(5/4));

            % Constant
            C_s_wale_2_num(ii,jj,kk) = sqrt(2).*S_ij_2.^(3/2);
            C_s_wale_2_den(ii,jj,kk) = S_ij_2.*S_ij_d_2.^(3/2)./(S_ij_2.^(5/2) + S_ij_d_2.^(5/4));

        end
    end
end


% % Div tensor
% divA_filt_favre = A_uu_filt_favre.^2 + A_vv_filt_favre.^2 + A_ww_filt_favre.^2;
% 
% % Traceless symmetric part of the square of the velocity gradient tensor
% S_d_uu = (A_uu_filt_favre.^2 + A_uu_filt_favre.^2) - 2/3*divA_filt_favre;
% S_d_uv = (A_uv_filt_favre.^2 + A_vu_filt_favre.^2);
% S_d_uw = (A_uw_filt_favre.^2 + A_wu_filt_favre.^2);
% S_d_vv = (A_vv_filt_favre.^2 + A_vv_filt_favre.^2) - 2/3*divA_filt_favre;
% S_d_vw = (A_vw_filt_favre.^2 + A_wv_filt_favre.^2);
% S_d_ww = (A_ww_filt_favre.^2 + A_ww_filt_favre.^2) - 2/3*divA_filt_favre;
% 
% % Inner product
% S_ij_d_2 = (0.5*S_d_uu).^2 + 2*(0.5*S_d_uv).^2 + 2*(0.5*S_d_uw).^2 + (0.5*S_d_vv).^2 + 2*(0.5*S_d_vw).^2 + (0.5*S_d_ww).^2;
% S_ij_2   = (0.5*S_uu_filt_favre).^2 + 2*(0.5*S_uv_filt_favre).^2 + 2*(0.5*S_uw_filt_favre).^2 + (0.5*S_vv_filt_favre).^2 + 2*(0.5*S_vw_filt_favre).^2 + (0.5*S_ww_filt_favre).^2;

% % Cw based on Smagorinsky (Cw ~ 0.55/0.6 for Cs ~ 0.18)
% C_s_wale_2_1 = Spatial_avg_XZ_var(sqrt(2).*S_ij_2.^(3/2),num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_2_2 = Spatial_avg_XZ_var(S_ij_2.*S_ij_d_2.^(3/2)./(S_ij_2.^(5/2) + S_ij_d_2.^(5/4)),num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_2   = C_s.^2.*C_s_wale_2_1./C_s_wale_2_2;

% C_w based on grid by grid point
C_s_wale_2_1 = Spatial_avg_XZ_var(C_s_wale_2_num,num_points_x,num_points_y,num_points_z, fi, idx);
C_s_wale_2_2 = Spatial_avg_XZ_var(C_s_wale_2_den,num_points_x,num_points_y,num_points_z, fi, idx);
C_s_wale_2   = C_s.^2.*C_s_wale_2_1./C_s_wale_2_2;

C_s_wale = 0;
y_norm   = -1*(y(2,idx_filt,2) - y(2,num_points_y-(idx_filt-1),2));

for jj = idx_filt:num_points_y-(idx_filt-1)
    delta_y = 0.5*(y(2,jj+1,2) - y(2,jj-1,2));
    C_s_wale  = C_s_wale  + sqrt(C_s_wale_2(jj))*delta_y/y_norm;
end


% % Eddy-viscosity
% ve_WALE = (C_s_wale.*delta_filt).^2.*(S_ij_d_2.^(3/2))./(S_ij_2.^(5/2) + S_ij_d_2.^(5/4));

% Eddy-viscosity
ve_WALE = (C_s_wale.*delta_filt).^2.*D_WALE;


% OP_1_uu      = (1/2*delta_uu.*1/2.*delta_uu).^(3/2);
% OP_2_uu      = ((1/2*S_uu_filt_favre*1/2.*S_uu_filt_favre).^(5/2) + (1/2*delta_uu.*1/2.*delta_uu).^(5/4));
% OP_1_uv      = (1/2*delta_uv.*1/2.*delta_uv).^(3/2);
% OP_2_uv      = ((1/2*S_uv_filt_favre.*1/2.*S_uv_filt_favre).^(5/2) + (1/2*delta_uv.*1/2.*delta_uv).^(5/4));
% OP_1_uw      = (1/2*delta_uw.*1/2.*delta_uw).^(3/2);
% OP_2_uw      = ((1/2*S_uw_filt_favre.*1/2.*S_uw_filt_favre).^(5/2) + (1/2*delta_uw.*1/2.*delta_uw).^(5/4));
% OP_1_vv      = (1/2*delta_vv.*1/2.*delta_vv).^(3/2);
% OP_2_vv      = ((1/2*S_vv_filt_favre.*1/2.*S_vv_filt_favre).^(5/2) + (1/2*delta_vv.*1/2.*delta_vv).^(5/4));
% OP_1_vw      = (1/2*delta_vw.*1/2.*delta_vw).^(3/2);
% OP_2_vw      = ((1/2*S_vw_filt_favre.*1/2.*S_vw_filt_favre).^(5/2) + (1/2*delta_vw.*1/2.*delta_vw).^(5/4));
% OP_1_ww      = (1/2*delta_ww.*1/2.*delta_ww).^(3/2);
% OP_2_ww      = ((1/2*S_ww_filt_favre.*1/2.*S_ww_filt_favre).^(5/2) + (1/2*delta_ww.*1/2.*delta_ww).^(5/4));

% C_s_wale_uu_2_1 = Spatial_avg_XZ_var(sqrt(2).*(1/2*S_uu_filt_favre.*1/2.*S_uu_filt_favre).^(3/2),num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_uu_2_2 = Spatial_avg_XZ_var(1/2*S_uu_filt_favre.*1/2.*S_uu_filt_favre.*OP_1_uu./OP_2_uu,num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_uu_2   = C_s.^2.*C_s_wale_uu_2_1./C_s_wale_uu_2_2;
% 
% C_s_wale_uv_2_1 = Spatial_avg_XZ_var(sqrt(2).*(1/2*S_uv_filt_favre.*1/2.*S_uv_filt_favre).^(3/2),num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_uv_2_2 = Spatial_avg_XZ_var(1/2*S_uv_filt_favre.*1/2.*S_uv_filt_favre.*OP_1_uv./OP_2_uv,num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_uv_2 = C_s.^2.*C_s_wale_uv_2_1./C_s_wale_uv_2_2;
% 
% C_s_wale_uw_2_1 = Spatial_avg_XZ_var(sqrt(2).*(1/2*S_uw_filt_favre.*1/2.*S_uw_filt_favre).^(3/2),num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_uw_2_2 = Spatial_avg_XZ_var(1/2*S_uw_filt_favre.*1/2.*S_uw_filt_favre.*OP_1_uw./OP_2_uw,num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_uw_2   = C_s.^2.*C_s_wale_uw_2_1./C_s_wale_uw_2_2;
% 
% C_s_wale_vv_2_1 = Spatial_avg_XZ_var(sqrt(2).*(1/2*S_vv_filt_favre.*1/2.*S_vv_filt_favre).^(3/2),num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_vv_2_2 = Spatial_avg_XZ_var(1/2*S_vv_filt_favre.*1/2.*S_vv_filt_favre.*OP_1_vv./OP_2_vv,num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_vv_2   = C_s.^2.*C_s_wale_vv_2_1./C_s_wale_vv_2_2;
% 
% C_s_wale_vw_2_1 = Spatial_avg_XZ_var(sqrt(2).*(1/2*S_vw_filt_favre.*1/2.*S_vw_filt_favre).^(3/2),num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_vw_2_2 = Spatial_avg_XZ_var(1/2*S_vw_filt_favre.*1/2.*S_vw_filt_favre.*OP_1_vw./OP_2_vw,num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_vw_2   = C_s.^2.*C_s_wale_vw_2_1./C_s_wale_vw_2_2;
% 
% C_s_wale_ww_2_1 = Spatial_avg_XZ_var(sqrt(2).*(1/2*S_ww_filt_favre.*1/2.*S_ww_filt_favre).^(3/2),num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_ww_2_2 = Spatial_avg_XZ_var(1/2*S_ww_filt_favre.*1/2.*S_ww_filt_favre.*OP_1_ww./OP_2_ww,num_points_x,num_points_y,num_points_z, fi, idx);
% C_s_wale_ww_2   = C_s.^2.*C_s_wale_ww_2_1./C_s_wale_ww_2_2;
% 
% C_s_wale_uu = 0; C_s_wale_uv = 0; C_s_wale_uw = 0; C_s_wale_vv = 0; C_s_wale_vw = 0; C_s_wale_ww = 0;
% 
% y_norm = y(2,idx_filt,2) - y(2,num_points_y-(idx_filt-1),2);
% 
% for jj = idx_filt:num_points_y-(idx_filt-1)
%     delta_y = 0.5*(y(2,jj+1,2) - y(2,jj-1,2));
%     C_s_wale_uu  = C_s_wale_uu  + C_s_wale_uu_2(jj)*delta_y/y_norm;
%     C_s_wale_uv  = C_s_wale_uv  + C_s_wale_uv_2(jj)*delta_y/y_norm;
%     C_s_wale_uw  = C_s_wale_uw  + C_s_wale_uw_2(jj)*delta_y/y_norm;
%     C_s_wale_vv  = C_s_wale_vv  + C_s_wale_vv_2(jj)*delta_y/y_norm;
%     C_s_wale_vw  = C_s_wale_vw  + C_s_wale_vw_2(jj)*delta_y/y_norm;
%     C_s_wale_ww  = C_s_wale_ww  + C_s_wale_ww_2(jj)*delta_y/y_norm;
% end
% 
% C_s_wale = sqrt(mean([C_s_wale_uu,C_s_wale_uv,C_s_wale_uw,C_s_wale_uv,C_s_wale_vv,C_s_wale_vw,C_s_wale_uw,C_s_wale_vw,C_s_wale_ww]));

% ve_WALE_uu = (C_s_wale.*delta_filt).^2.*OP_1_uu./OP_2_uu;
% ve_WALE_uv = (C_s_wale.*delta_filt).^2.*OP_1_uv./OP_2_uv;
% ve_WALE_uw = (C_s_wale.*delta_filt).^2.*OP_1_uw./OP_2_uw;
% ve_WALE_vv = (C_s_wale.*delta_filt).^2.*OP_1_vv./OP_2_vv;
% ve_WALE_vw = (C_s_wale.*delta_filt).^2.*OP_1_vw./OP_2_vw;
% ve_WALE_ww = (C_s_wale.*delta_filt).^2.*OP_1_ww./OP_2_ww;


Tau_xx_SFS{4} = - ve_WALE.*S_uu_filt_favre;
Tau_xy_SFS{4} = - ve_WALE.*S_uv_filt_favre;
Tau_xz_SFS{4} = - ve_WALE.*S_uw_filt_favre;
Tau_yy_SFS{4} = - ve_WALE.*S_vv_filt_favre;
Tau_yz_SFS{4} = - ve_WALE.*S_vw_filt_favre;
Tau_zz_SFS{4} = - ve_WALE.*S_ww_filt_favre;

%% Sigma Model (WALE Validation ONLY)
% Sigma constant 1.5 from Toda CTR Brief
C_Sigma = 1.5;

% Initialize D_Sigma operator
D_Sigma = zeros(size(dx));

% Compute for each gridpoint
for jj = idx_filt:num_points_x-(idx_filt-1)
    for ii = idx_filt:num_points_x-(idx_filt-1) % Only filtered points from y_outer
        for kk = idx_filt:num_points_z-(idx_filt-1)

            A_ij = [A_uu_filt_favre(ii,jj,kk) A_uv_filt_favre(ii,jj,kk) A_uw_filt_favre(ii,jj,kk);
                    A_uv_filt_favre(ii,jj,kk) A_vv_filt_favre(ii,jj,kk) A_vw_filt_favre(ii,jj,kk);
                    A_uw_filt_favre(ii,jj,kk) A_vw_filt_favre(ii,jj,kk) A_ww_filt_favre(ii,jj,kk)];

            A_ij2 = A_ij'*A_ij; % G = g^T g

            % Eigen values velocity gradient squared matrix (+ve)
            [A_2_eigenvec_raw, A_2_D_eigenval]   = eig(A_ij2);
             A_2_eigenval_raw                    = eig(A_ij2);
             % Squared root
             A_eigenval_raw                      = sqrt(A_2_eigenval_raw);
            [A_eigenval, idx_sort]               = sort(A_eigenval_raw,'descend'); % First greatest eigenvalue
             A_2_eigenvec                        = A_2_eigenvec_raw(:, idx_sort);    % Sort based on eigenvalue
             
            % D_Sigma
            D_Sigma(ii,jj,kk) = A_eigenval(3).*(A_eigenval(1) - A_eigenval(2)).*(A_eigenval(2) - A_eigenval(3))./(A_eigenval(1).^2);

        end
    end
end

% Eddy-viscosity model
ve_Sigma = (C_Sigma.*delta_filt).^2.*D_Sigma;

Tau_xx_SFS{6} = - ve_Sigma.*S_uu_filt_favre;
Tau_xy_SFS{6} = - ve_Sigma.*S_uv_filt_favre;
Tau_xz_SFS{6} = - ve_Sigma.*S_uw_filt_favre;
Tau_yy_SFS{6} = - ve_Sigma.*S_vv_filt_favre;
Tau_yz_SFS{6} = - ve_Sigma.*S_vw_filt_favre;
Tau_zz_SFS{6} = - ve_Sigma.*S_ww_filt_favre;


%% Similarity model
rhouu_filt_filt_favre    = FilterFields(rho_filt.*u_filt_favre.*u_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rhouv_filt_filt_favre    = FilterFields(rho_filt.*u_filt_favre.*v_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rhouw_filt_filt_favre    = FilterFields(rho_filt.*u_filt_favre.*w_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rhovv_filt_filt_favre    = FilterFields(rho_filt.*v_filt_favre.*v_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rhovw_filt_filt_favre    = FilterFields(rho_filt.*v_filt_favre.*w_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rhoww_filt_filt_favre    = FilterFields(rho_filt.*w_filt_favre.*w_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);

rhou_filt_filt_favre    = FilterFields(rho_filt.*u_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rhov_filt_filt_favre    = FilterFields(rho_filt.*v_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rhow_filt_filt_favre    = FilterFields(rho_filt.*w_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);

rho_filt_filt           = FilterFields(rho_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);

Tau_xx_SFS{5} = (rhouu_filt_filt_favre - rhou_filt_filt_favre.*rhou_filt_filt_favre./rho_filt_filt)./rho_filt;
Tau_xy_SFS{5} = (rhouv_filt_filt_favre - rhou_filt_filt_favre.*rhov_filt_filt_favre./rho_filt_filt)./rho_filt;
Tau_xz_SFS{5} = (rhouw_filt_filt_favre - rhou_filt_filt_favre.*rhow_filt_filt_favre./rho_filt_filt)./rho_filt;
Tau_yy_SFS{5} = (rhovv_filt_filt_favre - rhov_filt_filt_favre.*rhov_filt_filt_favre./rho_filt_filt)./rho_filt;
Tau_yz_SFS{5} = (rhovw_filt_filt_favre - rhov_filt_filt_favre.*rhow_filt_filt_favre./rho_filt_filt)./rho_filt;
Tau_zz_SFS{5} = (rhoww_filt_filt_favre - rhow_filt_filt_favre.*rhow_filt_filt_favre./rho_filt_filt)./rho_filt;


%% Models comparison
figure
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy_SFS{1}(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy_SFS{2}(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0.8500 0.3250 0.0980])
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy_SFS{3}(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0.3010 0.7450 0.9330])
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy_SFS{4}(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0.4940 0.1840 0.5560])
plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy_SFS{5}(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0.6350 0.0780 0.1840])
xlim([0 2])
legend('DNS filt','Smagorisnky','Dynami','AMD','WALE','Similarity')
xlabel('${y/\delta}$','interpreter','latex')
ylabel('${\tau_{yy}}$','interpreter','latex')
pbaspect([1.8 1 1])
legend('Location','northwest','box','off','NumColumns', 3)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]);

%% Correlation factors

% Y - positions
y_outer_vec  = [0.2 1 1.8]*delta_h;
idx_lim      = idx_filt:num_points_x-(idx_filt-1); % Only outer points for filt
C_corr       = zeros([length(y_outer_vec),length(Tau_xx_SFS)]);
C_corr_trace = zeros([length(y_outer_vec),length(Tau_kk_SFS)]);

for yy = 1:length(y_outer_vec)

    [~,idx_y]    = min(abs(y(1,:,1)-y_outer_vec(yy)));

    % Decomposition point by point
    for nModel = 1:length(Tau_xx_SFS)

        C_corr_uu  = sum(dot(Tau_xx_DNS_filt(idx_lim,idx_y,idx_lim),Tau_xx_SFS{nModel}(idx_lim,idx_y,idx_lim)))./(sum(dot(Tau_xx_DNS_filt(idx_lim,idx_y,idx_lim),Tau_xx_DNS_filt(idx_lim,idx_y,idx_lim))))^0.5./(sum(dot(Tau_xx_SFS{nModel}(idx_lim,idx_y,idx_lim),Tau_xx_SFS{nModel}(idx_lim,idx_y,idx_lim))))^0.5;
        C_corr_uv  = sum(dot(Tau_xy_DNS_filt(idx_lim,idx_y,idx_lim),Tau_xy_SFS{nModel}(idx_lim,idx_y,idx_lim)))./(sum(dot(Tau_xy_DNS_filt(idx_lim,idx_y,idx_lim),Tau_xy_DNS_filt(idx_lim,idx_y,idx_lim))))^0.5./(sum(dot(Tau_xy_SFS{nModel}(idx_lim,idx_y,idx_lim),Tau_xy_SFS{nModel}(idx_lim,idx_y,idx_lim))))^0.5;
        C_corr_uw  = sum(dot(Tau_xz_DNS_filt(idx_lim,idx_y,idx_lim),Tau_xz_SFS{nModel}(idx_lim,idx_y,idx_lim)))./(sum(dot(Tau_xz_DNS_filt(idx_lim,idx_y,idx_lim),Tau_xz_DNS_filt(idx_lim,idx_y,idx_lim))))^0.5./(sum(dot(Tau_xz_SFS{nModel}(idx_lim,idx_y,idx_lim),Tau_xz_SFS{nModel}(idx_lim,idx_y,idx_lim))))^0.5;
        C_corr_vv  = sum(dot(Tau_yy_DNS_filt(idx_lim,idx_y,idx_lim),Tau_yy_SFS{nModel}(idx_lim,idx_y,idx_lim)))./(sum(dot(Tau_yy_DNS_filt(idx_lim,idx_y,idx_lim),Tau_yy_DNS_filt(idx_lim,idx_y,idx_lim))))^0.5./(sum(dot(Tau_yy_SFS{nModel}(idx_lim,idx_y,idx_lim),Tau_yy_SFS{nModel}(idx_lim,idx_y,idx_lim))))^0.5;
        C_corr_vw  = sum(dot(Tau_yz_DNS_filt(idx_lim,idx_y,idx_lim),Tau_yz_SFS{nModel}(idx_lim,idx_y,idx_lim)))./(sum(dot(Tau_yz_DNS_filt(idx_lim,idx_y,idx_lim),Tau_yz_DNS_filt(idx_lim,idx_y,idx_lim))))^0.5./(sum(dot(Tau_yz_SFS{nModel}(idx_lim,idx_y,idx_lim),Tau_yz_SFS{nModel}(idx_lim,idx_y,idx_lim))))^0.5;
        C_corr_ww  = sum(dot(Tau_zz_DNS_filt(idx_lim,idx_y,idx_lim),Tau_zz_SFS{nModel}(idx_lim,idx_y,idx_lim)))./(sum(dot(Tau_zz_DNS_filt(idx_lim,idx_y,idx_lim),Tau_zz_DNS_filt(idx_lim,idx_y,idx_lim))))^0.5./(sum(dot(Tau_zz_SFS{nModel}(idx_lim,idx_y,idx_lim),Tau_zz_SFS{nModel}(idx_lim,idx_y,idx_lim))))^0.5;

        C_corr(yy,nModel) = mean([C_corr_uu,C_corr_uv,C_corr_uw,C_corr_uv,C_corr_vv,C_corr_vw,C_corr_uw,C_corr_vw,C_corr_ww]);
    end

    % Decomposition point by point
    for nModel_kk = 1:length(Tau_kk_SFS)

        C_corr_kk  = sum(dot(Tau_kk_DNS_filt(idx_lim,idx_y,idx_lim),Tau_kk_SFS{nModel_kk}(idx_lim,idx_y,idx_lim)))./(sum(dot(Tau_kk_DNS_filt(idx_lim,idx_y,idx_lim),Tau_kk_DNS_filt(idx_lim,idx_y,idx_lim))))^0.5./(sum(dot(Tau_kk_SFS{nModel_kk}(idx_lim,idx_y,idx_lim),Tau_kk_SFS{nModel_kk}(idx_lim,idx_y,idx_lim))))^0.5;
 
        C_corr_trace(yy,nModel_kk) = C_corr_kk;
    end

end


%% PDF Trace
PDF_Trace_SFS_Model(x,y,z,idx_filt,Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz, Tau_kk_SFS,u_b,delta_h,fi)


%% Eigen decomposition comparison
SFS_model = {'Smagorinsky', 'Dynamic', 'AMD', 'WALE', 'Similarity'};
for nModel = 1:length(Tau_xx_SFS) % Normalize by Vreman trace model
    EigenDecomposition_SFS_Model(x,y,z,idx_filt, Tau_xx_SFS{nModel},Tau_xy_SFS{nModel},Tau_xz_SFS{nModel},Tau_yy_SFS{nModel},Tau_yz_SFS{nModel},Tau_zz_SFS{nModel},Tau_kk_u_SFS{2},Tau_kk_v_SFS{2},Tau_kk_w_SFS{2},Tau_kk_SFS{2},u_b,delta_h, fi, SFS_model{nModel});
end



%% Equation of State Model

%% Taylor expansion model
J_t2.dcv_dT = 0;
J_t         = Jacobian_thermodynamics(bSolver, rho_filt,T_filt,P_filt, J_t2,Substance, HP_model);
d2P_drhodT  = J_t.d2P_drhodT;

rhoT_filt   = FilterFields(rho_filt.*T_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
delta_P     = 0.5*d2P_drhodT.*(rho_filt.*T_filt - rhoT_filt);
P_Taylor    = P_LES + delta_P;

P_Taylor_avg = Spatial_avg_XZ_var(P_Taylor,num_points_x,num_points_y,num_points_z, fi, idx);


%% ILA model
P_LES_filt   = FilterFields(P_LES,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
T_filt_filt  = FilterFields(T_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);

P_filt_LES   = Calculate_P_PengRobinson(rho_filt_filt,T_filt_filt,Substance);


T_filt_FG               = FilterFields(T_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
T_filt_FG_filt          = FilterFields(T_filt_FG,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
T_filt_FG_filt_FG       = FilterFields(T_filt_FG_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
rho_filt_FG_filt        = FilterFields(rho_filt_FG,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rho_filt_FG_filt_FG     = FilterFields(rho_filt_FG_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);

P_filt_FG_LES = Calculate_P_PengRobinson(rho_filt_FG,T_filt_FG,Substance);
P_filt_FG_LES_filt      = FilterFields(P_filt_FG_LES,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
P_filt_FG_LES_filt_FG   = FilterFields(P_filt_FG_LES_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);

P_filt_FG_filt_FG_LES   = Calculate_P_PengRobinson(rho_filt_FG_filt_FG,T_filt_FG_filt_FG,Substance);


P_LES_filt_FG            = FilterFields(P_LES_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
P_filt_LES_FG            = FilterFields(P_filt_LES,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);

M_P = (P_filt_FG_LES_filt_FG - P_filt_FG_filt_FG_LES) - (P_LES_filt_FG - P_filt_LES_FG);


P_LES_FG   = FilterFields(P_LES,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);

L_P         = P_LES_FG - P_filt_FG_LES;

% Coefficient C_p 
C_P = Spatial_avg_XZ_var(L_P.*M_P,num_points_x,num_points_y,num_points_z, fi, idx_filt)./Spatial_avg_XZ_var(M_P.*M_P,num_points_x,num_points_y,num_points_z, fi, idx_filt);
C_P(C_P<0) = 0;


P_ILA = P_LES + C_P.*(P_LES_filt - P_filt_LES);
P_ILA_avg = Spatial_avg_XZ_var(P_ILA,num_points_x,num_points_y,num_points_z, fi, idx);

% ,[0.4660 0.6740 0.1880])
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy_SFS{2}(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0.8500 0.3250 0.0980])
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy_SFS{3}(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0.3010 0.7450 0.9330])
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy_SFS{4}(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0.4940 0.1840 0.5560])
% plot(y(2,idx_filt:end-(idx_filt),2)/delta_h, Tau_yy_SFS{5}(idx_filt,idx_filt:end-(idx_filt),idx_filt),'linewidth',2, 'LineStyle','-','color',[0.6350 0.0780 0.1840])

% EOS
figure; hold on; box on
% set(gca,'yscale','log')
semilogy(y(2,idx:end-idx-1,2)/delta_h,(P_LES_avg(idx:end-idx-1))*Norm{6},'linewidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
hold on
semilogy(y(2,idx:end-idx-1,2)/delta_h,(P_filt_avg(idx:end-idx-1))*Norm{6},'linewidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
semilogy(y(2,idx:end-idx-1,2)/delta_h,(P_Taylor_avg(idx:end-idx-1))*Norm{6},'linewidth',2, 'LineStyle','-.','color',[0.3010 0.7450 0.9330])
semilogy(y(2,idx:end-idx-1,2)/delta_h,(P_ILA_avg(idx:end-idx-1))*Norm{6},'linewidth',2, 'LineStyle','--','color',[0.8500 0.3250 0.0980])

% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_6_avg{1}(idx:end-idx-1),'linewidth',2)
xlabel('${y/\delta}$','interpreter','latex')
ylabel('${\langle P^\ast \rangle}$','interpreter','latex')
legend([{'$P(\overline{\rho},\breve{T})^\ast$'},{'$\overline{P(\rho,T)}^\ast$'},{'${P(\overline{\rho},\breve{T})}^\ast + {\delta_P}^\ast$'},{'${P_{ILA}}^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.8 1 1])
legend('Location','southeast','box','off','NumColumns', 3)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]); %ylim([5900 6500])
% axis tight
% saveas(gca,strcat('Figures/Unclosed_EOS_', num2str(fi), 'xDelta_Comparison'),'png')
% saveas(gca,strcat('Figures/Unclosed_EOS_', num2str(fi), 'xDelta_Comparison'),'epsc')
exportgraphics(gca,strcat('Figures/Unclosed_EOS_', num2str(fi), 'xDelta_Comparison', '.jpeg'),'Resolution',300)
exportgraphics(gca,strcat('Figures/Unclosed_EOS_', num2str(fi), 'xDelta_Comparison', '.png'),'Resolution',300)


%% Correlation factors

% Y - positions
y_outer_vec  = [0.2 1 1.8]*delta_h;
C_corr_Taylor       = zeros(1,length(y_outer_vec));
C_corr_ILA          = zeros(1,length(y_outer_vec));

for yy = 1:length(y_outer_vec)

    [~,idx_y]    = min(abs(y(1,:,1)-y_outer_vec(yy)));


    C_corr_Taylor(yy)  = sum(dot(P_filt(idx_lim,idx_y,idx_lim),P_Taylor(idx_lim,idx_y,idx_lim)))./(sum(dot(P_filt(idx_lim,idx_y,idx_lim),P_filt(idx_lim,idx_y,idx_lim))))^0.5./(sum(dot(P_Taylor(idx_lim,idx_y,idx_lim),P_Taylor(idx_lim,idx_y,idx_lim))))^0.5;
    C_corr_ILA(yy)     = sum(dot(P_filt(idx_lim,idx_y,idx_lim),P_ILA(idx_lim,idx_y,idx_lim)))./(sum(dot(P_filt(idx_lim,idx_y,idx_lim),P_filt(idx_lim,idx_y,idx_lim))))^0.5./(sum(dot(P_ILA(idx_lim,idx_y,idx_lim),P_ILA(idx_lim,idx_y,idx_lim))))^0.5;


end
