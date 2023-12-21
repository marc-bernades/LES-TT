%% Build the unclosed terms based on DNS
clear all; close all; clc

%% Load filtered data (ensemble-average)
% Name of the .mat file containing the filtered ensemble-average data for
% filter widths to analyse
name_file_out_filt_dataset = 'Filtered_DNS_ensemble_average';
% Load mat file
Load_EnsembleAverage_Filt;

%% Problem fixed parameters
delta_h   = 100*1E-6; % Channel half height
bSolver           = 'Real';                    % Define Ideal or Real Gas Model
Substance         = 'N2';                      % Fluid selected substance
HP_model          = 'HighPressure';            % Transport coefficients model: 'Constant', 'LowPressure', 'HighPressure'
Fluid.R_specific  = 287.058;
bPlot = false; % Plot boolean
fi    = 8;     % Filter width chosen for the results

%% Obtain tau and bulk metrics from average DNS 
% Select file
source_path  = '/home/marc/Documents/Doctorat/Solver_and_postprocess/RHEA/Paper3/Post_process/';
File_name    = '3d_high_pressure_turbulent_channel_flow_69100000.h5';
File_name_complete = strcat(source_path,File_name);
Info      = h5info(File_name_complete);
Name_var  = {Info.Datasets.Name};

% Load variables
for ii = 1:length(Name_var)
    value = h5read(File_name_complete,strcat('/',Name_var{ii}) );
    assignin('base',Name_var{ii}, value)
    %     x = h5read('3d_high_pressure_turbulent_channel_flow_58900000.h5','/x');
end

% Add fields
nu       = mu./rho;
avg_nu   = avg_mu./avg_rho;

% Grid points
num_points_x  = length(x(:,1,1));
num_points_y  = length(x(1,:,1));
num_points_z  = length(x(1,1,:));
num_points_xz = num_points_x*num_points_z;

% Delta x, y and z based on CentralDerivative_d1_2ndOrder
[~,dx,~] = CentralDerivative_d1_2ndOrder(x);
[dy,~,~] = CentralDerivative_d1_2ndOrder(y);
[~,~,dz] = CentralDerivative_d1_2ndOrder(z);

% Bulk and tau values
[rho_b, mu_b, P_b, T_b, u_b, Re_b, volume, delta] = Calculate_bulk_delta_metrics(x,y,z,rho, mu, P, T, u, num_points_x, num_points_y, num_points_z, delta_h, dx, dy, dz);
[rho_bw, mu_bw, T_tau_bw, T_bw, u_tau_bw,rho_tw, mu_tw, T_tau_tw, T_tw, u_tau_tw] = Calculate_wall_metrics(x,y,z,avg_rho,avg_mu,avg_T,avg_c_p,avg_kappa,avg_u,rho_b,T_b,u_b,num_points_x,num_points_y,num_points_z,delta_h);
[y_plus_bw,u_plus_bw,T_plus_bw,y_plus_tw,u_plus_tw,T_plus_tw] = Transform_WallUnits(y,avg_u,avg_T,u_tau_bw,rho_bw,mu_bw,T_bw,T_tau_bw,u_tau_tw,rho_tw,mu_tw,T_tw,T_tau_tw,num_points_x,num_points_y,num_points_z,delta_h);



%% Create fields
rho_filt       = eval(strcat('rho','_filt',num2str(fi),'xDelta'));
u_filt         = eval(strcat('u','_filt',num2str(fi),'xDelta'));
v_filt         = eval(strcat('v','_filt',num2str(fi),'xDelta'));
w_filt         = eval(strcat('w','_filt',num2str(fi),'xDelta'));
P_filt         = eval(strcat('P','_filt',num2str(fi),'xDelta'));
T_filt         = eval(strcat('T','_filt',num2str(fi),'xDelta'));
mu_filt        = eval(strcat('mu','_filt',num2str(fi),'xDelta'));
kappa_filt     = eval(strcat('kappa','_filt',num2str(fi),'xDelta'));
c_p_filt       = eval(strcat('c_p','_filt',num2str(fi),'xDelta'));
c_v_filt       = eval(strcat('c_v','_filt',num2str(fi),'xDelta'));
u_filt_favre   = eval(strcat('u_filt_favre',num2str(fi),'xDelta'));
uu_filt_favre  = eval(strcat('uu_filt_favre',num2str(fi),'xDelta'));
v_filt_favre   = eval(strcat('v_filt_favre',num2str(fi),'xDelta'));
vv_filt_favre  = eval(strcat('vv_filt_favre',num2str(fi),'xDelta'));
w_filt_favre   = eval(strcat('w_filt_favre',num2str(fi),'xDelta'));
ww_filt_favre  = eval(strcat('ww_filt_favre',num2str(fi),'xDelta'));
uv_filt_favre  = eval(strcat('uv_filt_favre',num2str(fi),'xDelta'));
uw_filt_favre  = eval(strcat('uw_filt_favre',num2str(fi),'xDelta'));
vw_filt_favre  = eval(strcat('vw_filt_favre',num2str(fi),'xDelta'));
rhou_visc_filt = eval(strcat('rhou_visc_filt',num2str(fi),'xDelta'));
rhov_visc_filt = eval(strcat('rhov_visc_filt',num2str(fi),'xDelta'));
rhow_visc_filt = eval(strcat('rhow_visc_filt',num2str(fi),'xDelta'));
u_grad_P_filt  = eval(strcat('u_grad_P_filt',num2str(fi),'xDelta'));
alpha_4_filt   = eval(strcat('alpha_4_filt',num2str(fi),'xDelta'));
alpha_5_filt   = eval(strcat('alpha_5_filt',num2str(fi),'xDelta'));

% avg_rho_xz_filt   = eval(strcat('avg_rho_xz_filt',num2str(fi),'xDelta'));
% avg_mu_xz_filt    = eval(strcat('avg_mu_xz_filt',num2str(fi),'xDelta'));
% avg_u_xz_filt     = eval(strcat('avg_u_xz_filt',num2str(fi),'xDelta'));
% avg_v_xz_filt     = eval(strcat('avg_v_xz_filt',num2str(fi),'xDelta'));
% avg_w_xz_filt     = eval(strcat('avg_w_xz_filt',num2str(fi),'xDelta'));
% avg_rhou_xz_filt  = eval(strcat('avg_rhou_xz_filt',num2str(fi),'xDelta'));
% avg_rhov_xz_filt  = eval(strcat('avg_rhov_xz_filt',num2str(fi),'xDelta'));
% avg_rhow_xz_filt  = eval(strcat('avg_rhow_xz_filt',num2str(fi),'xDelta'));



%% Filtered fields

% Permute matrix to have X on second position as HPCFS
rho_perm = permute(rho_filt, [2,1,3]);
u_perm   = permute(u_filt_favre,   [2,1,3]);
v_perm   = permute(v_filt_favre,   [2,1,3]);
w_perm   = permute(w_filt_favre,   [2,1,3]);

% Inviscid components
rho_conv  = Inviscid_CD(rho_perm,u_perm,v_perm,w_perm,ones(size(x)),dx,dy,dz);
rhou_conv = Inviscid_CD(rho_perm,u_perm,v_perm,w_perm,u_perm,dx,dy,dz);
rhov_conv = Inviscid_CD(rho_perm,u_perm,v_perm,w_perm,v_perm,dx,dy,dz);
rhow_conv = Inviscid_CD(rho_perm,u_perm,v_perm,w_perm,w_perm,dx,dy,dz);

% Permute back
rho_conv_filt  = permute(rho_conv,  [2,1,3]);
rhou_conv_filt = permute(rhou_conv, [2,1,3]);
rhov_conv_filt = permute(rhov_conv, [2,1,3]);
rhow_conv_filt = permute(rhow_conv, [2,1,3]);

%% Momentum
%% alpha 1
% Turbulent stress tensor
Tau_xx = (uu_filt_favre  - u_filt_favre.^2);
Tau_xy = (uv_filt_favre  - u_filt_favre.*v_filt_favre);
Tau_xz = (uw_filt_favre  - u_filt_favre.*w_filt_favre);
Tau_yy = (vv_filt_favre  - v_filt_favre.^2);
Tau_yz = (vw_filt_favre  - v_filt_favre.*w_filt_favre);
Tau_zz = (ww_filt_favre  - w_filt_favre.^2);

[~,d_alpha_xx_x,~] = CentralDerivative_d1_2ndOrder(rho_filt.*Tau_xx);
[d_alpha_yy_y,~,~] = CentralDerivative_d1_2ndOrder(rho_filt.*Tau_yy);
[~,~,d_alpha_zz_z] = CentralDerivative_d1_2ndOrder(rho_filt.*Tau_zz);

[d_alpha_xy_y,d_alpha_xy_x,~]            = CentralDerivative_d1_2ndOrder(rho_filt.*Tau_xy);
[d_alpha_xz_y,d_alpha_xz_x,d_alpha_xz_z] = CentralDerivative_d1_2ndOrder(rho_filt.*Tau_xz);
[d_alpha_yz_y,~,           d_alpha_yz_z] = CentralDerivative_d1_2ndOrder(rho_filt.*Tau_yz);


alpha_1{1} = d_alpha_xx_x./dx + d_alpha_xy_y./dy + d_alpha_xz_z./dz;
alpha_1{2} = d_alpha_xy_x./dx + d_alpha_yy_y./dy + d_alpha_yz_z./dz;
alpha_1{3} = d_alpha_xz_x./dx + d_alpha_yz_y./dy + d_alpha_zz_z./dz;


%% alpha_2 DNS: stress tensor (sigma)
[rhou_visc_LES, rhov_visc_LES, rhow_visc_LES, rhoE_visc_LES, Tau_LES, q_div_LES] = Calculate_Viscous_LES(u_filt_favre,v_filt_favre,w_filt_favre,mu_filt,kappa_filt,T_filt,dx,dy,dz,1,1);

alpha_2{1}    = rhou_visc_filt - rhou_visc_LES;
alpha_2{2}    = rhov_visc_filt - rhov_visc_LES;
alpha_2{3}    = rhow_visc_filt - rhow_visc_LES;

%% Pressure
%% alpha_3
[dP_y_filt,dP_x_filt,dP_z_filt]        = CentralDerivative_d1_2ndOrder(P_filt);
u_grad_P_LES     = u_filt_favre.*dP_x_filt./dx + v_filt_favre.*dP_y_filt./dy + w_filt_favre.*dP_z_filt./dz;

alpha_3{1}       = u_grad_P_LES - u_grad_P_filt;

%% alpha 4
[du_y_filt_favre,du_x_filt_favre,du_z_filt_favre]  = CentralDerivative_d1_2ndOrder(u_filt_favre);
[dv_y_filt_favre,dv_x_filt_favre,dv_z_filt_favre]  = CentralDerivative_d1_2ndOrder(v_filt_favre);
[dw_y_filt_favre,dw_x_filt_favre,dw_z_filt_favre]  = CentralDerivative_d1_2ndOrder(w_filt_favre);

div_U_filt_favre = du_x_filt_favre./dx + dv_y_filt_favre./dy + dw_z_filt_favre./dz;


[sos_filt, Beta_T_filt, Beta_v_filt, Beta_s_filt, Alpha_p_filt] = Calculate_sos(bSolver,rho_filt,T_filt,c_p_filt,P_filt,0,Fluid,Substance);

alpha_4_LES    = rho_filt.*sos_filt.^2.*div_U_filt_favre;

alpha_4{1}     = alpha_4_LES - alpha_4_filt;


%% alpha 5
Tau_dU_LES     = (du_x_filt_favre./dx).*Tau_LES.Tau_xx + (du_y_filt_favre./dy).*Tau_LES.Tau_xy + (du_z_filt_favre./dz).*Tau_LES.Tau_xz  + ...
               + (dv_x_filt_favre./dx).*Tau_LES.Tau_yx + (dv_y_filt_favre./dy).*Tau_LES.Tau_yy + (dv_z_filt_favre./dz).*Tau_LES.Tau_yz  + ...
               + (dw_x_filt_favre./dx).*Tau_LES.Tau_zx + (dw_y_filt_favre./dy).*Tau_LES.Tau_zy + (dw_z_filt_favre./dz).*Tau_LES.Tau_zz;

alpha_5_LES   = Alpha_p_filt./(c_v_filt.*Beta_T_filt).*(1./rho_filt.*(Tau_dU_LES - q_div_LES));

alpha_5{1}    = alpha_5_filt - alpha_5_LES;

%% Equation of state
%% alpha 6
P_LES        = Calculate_P_PengRobinson(rho_filt,T_filt,Substance);
alpha_6{1}   = P_filt - P_LES;



%% Ensembled average unclosed terms
Re_tau_bw = 100;
l_v       = mu_bw/rho_bw/u_tau_bw; % Viscous length scale
delta_BL  = Re_tau_bw*l_v;
y_delta_BL_target = 0.2; % Outer wall target (Larsson 2015)
y_outer           = y_delta_BL_target*delta_BL; 

[~,idx]    = min(abs(y(1,:,1)-y_outer));
idx_filt   = idx;
y_outer_norm = y(1,idx,1)/delta_h;
disp("Outer position at y/delta = " + num2str(y_outer_norm) + ...
    " and y_plus_bw = " + num2str(y_plus_bw(idx)) + " and y_plus_tw = " + num2str(y_plus_tw(idx)))



[alpha_1_avg{1}] = Spatial_avg_XZ_var(alpha_1{1},num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_1_avg{2}] = Spatial_avg_XZ_var(alpha_1{2},num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_1_avg{3}] = Spatial_avg_XZ_var(alpha_1{3},num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_2_avg{1}] = Spatial_avg_XZ_var(alpha_2{1},num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_2_avg{2}] = Spatial_avg_XZ_var(alpha_2{2},num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_2_avg{3}] = Spatial_avg_XZ_var(alpha_2{3},num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_3_avg{1}] = Spatial_avg_XZ_var(alpha_3{1},num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_4_avg{1}] = Spatial_avg_XZ_var(alpha_4{1},num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_5_avg{1}] = Spatial_avg_XZ_var(alpha_5{1},num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_6_avg{1}] = Spatial_avg_XZ_var(alpha_6{1},num_points_x,num_points_y,num_points_z, fi, idx);

Norm{1} = delta_h/(u_b^2*rho_b);
Norm{2} = delta_h^2/(u_b*mu_b);
Norm{3} = delta_h/(u_b^2*rho_b);
Norm{4} = Norm{3};
Norm{5} = Norm{3};
Norm{6} = 1/(u_b*rho_b);




%% Ensembled average - Momentum
% Convective
[rhou_conv_LES_avg] = Spatial_avg_XZ_var(rhou_conv_filt,num_points_x,num_points_y,num_points_z, fi, idx);
[rhov_conv_LES_avg] = Spatial_avg_XZ_var(rhov_conv_filt,num_points_x,num_points_y,num_points_z, fi, idx);
[rhow_conv_LES_avg] = Spatial_avg_XZ_var(rhow_conv_filt,num_points_x,num_points_y,num_points_z, fi, idx);
% Viscous
[rhou_visc_LES_avg] = Spatial_avg_XZ_var(rhou_visc_LES,num_points_x,num_points_y,num_points_z, fi, idx);
[rhov_visc_LES_avg] = Spatial_avg_XZ_var(rhov_visc_LES,num_points_x,num_points_y,num_points_z, fi, idx);
[rhow_visc_LES_avg] = Spatial_avg_XZ_var(rhow_visc_LES,num_points_x,num_points_y,num_points_z, fi, idx);
% Pressure
dP_rhou_LES = dP_x_filt./(dx);
dP_rhov_LES = dP_y_filt./(dy);
dP_rhow_LES = dP_z_filt./(dz);
[dP_rhou_LES_avg]    = Spatial_avg_XZ_var(dP_rhou_LES,num_points_x,num_points_y,num_points_z, fi, idx);
[dP_rhov_LES_avg]    = Spatial_avg_XZ_var(dP_rhov_LES,num_points_x,num_points_y,num_points_z, fi, idx);
[dP_rhow_LES_avg]    = Spatial_avg_XZ_var(dP_rhow_LES,num_points_x,num_points_y,num_points_z, fi, idx);

% Transient
rhou_trans = -rhou_conv_filt + rhou_visc_LES - dP_rhou_LES_avg - alpha_1{1} + alpha_2{1};
rhov_trans = -rhov_conv_filt + rhov_visc_LES - dP_rhov_LES_avg - alpha_1{2} + alpha_2{2};
rhow_trans = -rhow_conv_filt + rhow_visc_LES - dP_rhow_LES_avg - alpha_1{3} + alpha_2{3};
[rhou_trans_avg]    = Spatial_avg_XZ_var(rhou_trans,num_points_x,num_points_y,num_points_z, fi, idx);
[rhov_trans_avg]    = Spatial_avg_XZ_var(rhov_trans,num_points_x,num_points_y,num_points_z, fi, idx);
[rhow_trans_avg]    = Spatial_avg_XZ_var(rhow_trans,num_points_x,num_points_y,num_points_z, fi, idx);


%% Ensembled average - Pressure
[u_grad_P_LES_avg]  = Spatial_avg_XZ_var(u_grad_P_LES,num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_4_LES_avg]   = Spatial_avg_XZ_var(alpha_4_LES,num_points_x,num_points_y,num_points_z, fi, idx);
[alpha_5_LES_avg]   = Spatial_avg_XZ_var(alpha_5_LES,num_points_x,num_points_y,num_points_z, fi, idx);

% Transient
P_trans         = - u_grad_P_LES - alpha_4_LES + alpha_5_LES + alpha_3{1} + alpha_4{1} + alpha_5{1};
[P_trans_avg]   = Spatial_avg_XZ_var(P_trans,num_points_x,num_points_y,num_points_z, fi, idx);

%% Ensembled average - EOS
[P_LES_avg]   = Spatial_avg_XZ_var(P_LES,num_points_x,num_points_y,num_points_z, fi, idx);
[P_filt_avg]  = Spatial_avg_XZ_var(P_filt,num_points_x,num_points_y,num_points_z, fi, idx);


%% Pressure gradient volumetric force
rho_bw_filt = 0; rho_tw_filt = 0;
mu_bw_filt  = 0; mu_tw_filt  = 0;
u_inner_bw_filt    = 0; u_inner_tw_filt    = 0;
u_boundary_bw_filt = 0; u_boundary_tw_filt = 0;
a_t_bw = 0; a_t_tw = 0;

for ii = 2:num_points_x-1
    for kk = 2:num_points_z-1
        % Bottom wall
        jj = 1;
        a_bw = dx(ii,jj,kk)*dz(ii,jj,kk);
        a_t_bw = a_t_bw + a_bw;
        rho_bw_filt = rho_bw_filt + a_bw*0.5*(rho_filt(ii,jj,kk) + rho_filt(ii,jj+1,kk));
        mu_bw_filt  = mu_bw_filt  + a_bw*0.5*(mu_filt(ii,jj,kk)  + mu_filt(ii,jj+1,kk));
        u_boundary_bw_filt = u_boundary_bw_filt + a_bw*u_filt(ii,jj,kk);
        u_inner_bw_filt    = u_inner_bw_filt    + a_bw*u_filt(ii,jj+1,kk);

        % Top wall
        jj = num_points_y;
        a_tw = dx(ii,jj,kk)*dz(ii,jj,kk);
        a_t_tw = a_t_tw + a_tw;
        rho_tw_filt = rho_tw_filt + a_tw*0.5*(rho_filt(ii,jj,kk) + rho_filt(ii,jj-1,kk));
        mu_tw_filt  = mu_tw_filt  + a_tw*0.5*(mu_filt(ii,jj,kk)  + mu_filt(ii,jj-1,kk));
        u_boundary_tw_filt = u_boundary_tw_filt + a_tw*u_filt(ii,jj,kk);
        u_inner_tw_filt    = u_inner_tw_filt    + a_tw*u_filt(ii,jj-1,kk);

    end
end


rho_bw_filt  = rho_bw_filt/a_t_bw; rho_tw_filt  = rho_tw_filt/a_t_tw;
mu_bw_filt   = mu_bw_filt/a_t_bw;  mu_tw_filt  = mu_tw_filt/a_t_tw;
u_inner_bw_filt = u_inner_bw_filt/a_t_bw; u_inner_tw_filt = u_inner_tw_filt/a_t_tw; 
u_boundary_bw_filt = u_boundary_bw_filt/a_t_bw; u_boundary_tw_filt = u_boundary_tw_filt/a_t_tw; 
tau_bw = mu_bw_filt*(u_inner_bw_filt - u_boundary_bw_filt)/(y(1,2,1) - y(1,1,1));
tau_tw = mu_tw_filt*(u_inner_tw_filt - u_boundary_tw_filt)/(y(1,end,1) - y(1,end-1,1));

dP_vol = - (tau_tw + tau_bw)/(2*delta_h);
% dP_rhou_LES_avg = dP_rhou_LES_avg +  dP_vol;


%% PLOTS unclosed terms evolution and breakdown
if bPlot
    Plot_UnclosedTerms;
end

