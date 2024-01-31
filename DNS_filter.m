%% Load h5
clear all; close all; clc
tic

%% Save and plot data
bSave = true;
bPlot = false;

%% Problem fixed parameters
bSolver           = 'Real';                    % Define Ideal or Real Gas Model
% Substance
Substance         = 'N2';                      % Fluid selected substance
HP_model          = 'HighPressure';                % Transport coefficients model: 'Constant', 'LowPressure', 'HighPressure'
Fluid.R_specific  = 287.058;
% Dimensions
delta_h   = 100*1E-6; % Channel half height
L_x = 4*pi*delta_h;   % Length in the x direction
L_y = 2*delta_h;      % Length in the y direction
L_z = 4/3*pi*delta_h; % Length in the z direction

% Filter width vector
fi           = [1 2 4 6 8 10];

% Select file where DNS datasets are stored
source_path  = '/home/marc/Documents/Doctorat/Solver_and_postprocess/RHEA/Paper3/Post_process/';
% High-pressure
% File_name{1} = '3d_high_pressure_turbulent_channel_flow_49800000.h5';
% File_name{2} = '3d_high_pressure_turbulent_channel_flow_52100000.h5';
% File_name{3} = '3d_high_pressure_turbulent_channel_flow_54200000.h5';
% File_name{4} = '3d_high_pressure_turbulent_channel_flow_56200000.h5';
% File_name{5} = '3d_high_pressure_turbulent_channel_flow_58900000.h5';
% File_name{6} = '3d_high_pressure_turbulent_channel_flow_67100000.h5';
File_name{1} = '3d_high_pressure_turbulent_channel_flow_69100000.h5';
% Low-pressure 
% File_name = 'restart_data_file_LP.h5'; % Laminar
% File_name{1} = '3d_turbulent_channel_flow_56800000_LP_isoT.h5'; % Turbulent

n_datasets = length(File_name);

% Pre-allocate total TKE vector
TKE_tot_filt     = zeros(length(fi),n_datasets);

for n_data = 1:n_datasets

    disp("Computing data set " + File_name{n_data} + "...")

    File_name_complete = strcat(source_path,File_name{n_data});

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

    % Save for ensembled average
    u_avg_dataset(:,:,:,n_data)     = u;
    v_avg_dataset(:,:,:,n_data)     = v;
    w_avg_dataset(:,:,:,n_data)     = w;
    rho_avg_dataset(:,:,:,n_data)   = rho;
    P_avg_dataset(:,:,:,n_data)     = P;
    T_avg_dataset(:,:,:,n_data)     = T;
    mu_avg_dataset(:,:,:,n_data)    = mu;
    kappa_avg_dataset(:,:,:,n_data) = kappa;
    c_p_avg_dataset(:,:,:,n_data)   = c_p;
    c_v_avg_dataset(:,:,:,n_data)   = c_v;

    % Grid points
    num_points_x  = length(x(:,1,1));
    num_points_y  = length(x(1,:,1));
    num_points_z  = length(x(1,1,:));
    num_points_xz = num_points_x*num_points_z;

    % Delta x, y and z based on CentralDerivative_d1_2ndOrder
    [~,dx,~] = CentralDerivative_d1_2ndOrder(x);
    [dy,~,~] = CentralDerivative_d1_2ndOrder(y);
    [~,~,dz] = CentralDerivative_d1_2ndOrder(z);


    %% Post-process obtain fields
    % Calculate bulk values, cell size and deltas
    [rho_b(n_data), mu_b(n_data), P_b(n_data), T_b(n_data), u_b(n_data), Re_b(n_data), volume, delta] = Calculate_bulk_delta_metrics(x,y,z,rho, mu, P, T, u, num_points_x, num_points_y, num_points_z, delta_h, dx, dy, dz);

    % Calculate Favre fluctuations
    % First calculate average on XZ avg_rhoui, avg_rho, avg_ui
    [avg_rho_xz(n_data,:),avg_mu_xz(n_data,:),avg_u_xz(n_data,:),avg_v_xz(n_data,:),avg_w_xz(n_data,:),avg_rhou_xz(n_data,:),avg_rhov_xz(n_data,:),avg_rhow_xz(n_data,:)] = Spatial_avg_XZ(rho,mu,u,v,w,num_points_x,num_points_y,num_points_z);

    % Calculate favre fluctuations with spatial average
    [R_favre_uu_fluct(:,:,:,n_data),R_favre_vv_fluct(:,:,:,n_data),R_favre_ww_fluct(:,:,:,n_data)] = Calculate_favre_fluctuations_XZ(u,v,w,avg_rho_xz(n_data,:),avg_rhou_xz(n_data,:),avg_rhov_xz(n_data,:),avg_rhow_xz(n_data,:),num_points_y);

    % TKE and average Favre
    [TKE(:,:,:,n_data), TKE_y(n_data,:),TKE_tw(n_data,:),TKE_bw(n_data,:),avg_R_favre_uu_bw(n_data,:), avg_R_favre_vv_bw(n_data,:), avg_R_favre_ww_bw(n_data,:), ...
        avg_R_favre_uu_tw(n_data,:), avg_R_favre_vv_tw(n_data,:), avg_R_favre_ww_tw(n_data,:)] = Calculate_TKE(R_favre_uu_fluct(:,:,:,n_data),R_favre_vv_fluct(:,:,:,n_data),R_favre_ww_fluct(:,:,:,n_data),num_points_x,num_points_y,num_points_z);

    % Calculate wall values (based on time-averaged last dataset)
    [rho_bw, mu_bw, T_tau_bw, T_bw, u_tau_bw,rho_tw, mu_tw, T_tau_tw, T_tw, u_tau_tw] = Calculate_wall_metrics(x,y,z,avg_rho,avg_mu,avg_T,avg_c_p,avg_kappa,avg_u,rho_b(n_data),T_b(n_data),u_b(n_data),num_points_x,num_points_y,num_points_z,delta_h);

    % Normalize wall-units (based on time-averaged last dataset)
    [y_plus_bw,u_plus_bw,T_plus_bw,y_plus_tw,u_plus_tw,T_plus_tw] = Transform_WallUnits(y,avg_u,avg_T,u_tau_bw,rho_bw,mu_bw,T_bw,T_tau_bw,u_tau_tw,rho_tw,mu_tw,T_tw,T_tau_tw,num_points_x,num_points_y,num_points_z,delta_h);

    % TKE for DNS
    TKE_tot_filt(1,n_data)  = sum(TKE_y(n_data,2:end-1).*dy(2,2:end-1,2))/(2*delta_h); % Assign DNS to first position

    %% DNS Filter

    for nn = 2:length(fi)

        disp("Computing filter operations for filter width = " + num2str(fi(nn)) + " x delta...")

        % Filter width
        delta_filt = fi(nn)*delta;

        % Filter position
        n_filter = nn - 1;

        % Fields to filter
        Var_filt = {'u','v','w','rho','P','T','mu','E','nu', 'kappa', 'c_p', 'c_v'...
            'avg_u','avg_v','avg_w','avg_rho','avg_P','avg_T','avg_mu','avg_kappa','avg_c_p','avg_E','avg_nu','avg_rhou','avg_rhov','avg_rhow'};

        for ii = 1:length(Var_filt)
            value = eval(Var_filt{ii});

            % Choose filter
            Filter_type = 'CDLF_Box'; % 'Top_hat', 'Top_hat_2_Forward', 'Digital_image', 'CDLF_Lap', 'CDLF_Box'
            value_filt  = FilterFields(value,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type);
            assignin('base',strcat(Var_filt{ii},'_filt',num2str(fi(nn)),'xDelta'), value_filt)

        end

        % Filtered quantities
        u_filt(:,:,:,n_data,n_filter)          = eval(strcat('u','_filt',num2str(fi(nn)),'xDelta'));
        v_filt(:,:,:,n_data,n_filter)          = eval(strcat('v','_filt',num2str(fi(nn)),'xDelta'));
        w_filt(:,:,:,n_data,n_filter)          = eval(strcat('w','_filt',num2str(fi(nn)),'xDelta'));
        rho_filt(:,:,:,n_data,n_filter)        = eval(strcat('rho','_filt',num2str(fi(nn)),'xDelta'));
        P_filt(:,:,:,n_data,n_filter)          = eval(strcat('P','_filt',num2str(fi(nn)),'xDelta'));
        T_filt(:,:,:,n_data,n_filter)          = eval(strcat('T','_filt',num2str(fi(nn)),'xDelta'));
        mu_filt(:,:,:,n_data,n_filter)         = eval(strcat('mu','_filt',num2str(fi(nn)),'xDelta'));
        kappa_filt(:,:,:,n_data,n_filter)      = eval(strcat('kappa','_filt',num2str(fi(nn)),'xDelta'));
        c_p_filt(:,:,:,n_data,n_filter)        = eval(strcat('c_p','_filt',num2str(fi(nn)),'xDelta'));
        c_v_filt(:,:,:,n_data,n_filter)        = eval(strcat('c_v','_filt',num2str(fi(nn)),'xDelta'));


        % Additional filtered quantities to asssign
        u_filt_favre(:,:,:,n_data,n_filter)   = FilterFields(rho.*u,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type)./rho_filt(:,:,:,n_data,n_filter); 
        uu_filt_favre(:,:,:,n_data,n_filter)  = FilterFields(rho.*u.*u,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type)./rho_filt(:,:,:,n_data,n_filter); 
        v_filt_favre(:,:,:,n_data,n_filter)   = FilterFields(rho.*v,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type)./rho_filt(:,:,:,n_data,n_filter); 
        vv_filt_favre(:,:,:,n_data,n_filter)  = FilterFields(rho.*v.*v,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type)./rho_filt(:,:,:,n_data,n_filter); 
        w_filt_favre(:,:,:,n_data,n_filter)   = FilterFields(rho.*w,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type)./rho_filt(:,:,:,n_data,n_filter);
        ww_filt_favre(:,:,:,n_data,n_filter)  = FilterFields(rho.*w.*w,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type)./rho_filt(:,:,:,n_data,n_filter); 
        uv_filt_favre(:,:,:,n_data,n_filter)  = FilterFields(rho.*u.*v,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type)./rho_filt(:,:,:,n_data,n_filter); 
        uw_filt_favre(:,:,:,n_data,n_filter)  = FilterFields(rho.*u.*w,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type)./rho_filt(:,:,:,n_data,n_filter);
        vw_filt_favre(:,:,:,n_data,n_filter)  = FilterFields(rho.*v.*w,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type)./rho_filt(:,:,:,n_data,n_filter);

        [rhou_visc, rhov_visc, rhow_visc, rhoE_visc, Tau, q_div]  = Calculate_Viscous_LES(u,v,w,mu,kappa,T,dx,dy,dz,1,1);
        rhou_visc_filt(:,:,:,n_data,n_filter) = FilterFields(rhou_visc,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type);
        rhov_visc_filt(:,:,:,n_data,n_filter) = FilterFields(rhov_visc,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type);
        rhow_visc_filt(:,:,:,n_data,n_filter) = FilterFields(rhow_visc,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type);

        [dP_y,     dP_x,     dP_z]             = CentralDerivative_d1_2ndOrder(P);
        u_grad_P                               = u.*dP_x./dx + v.*dP_y./dy + w.*dP_z./dz;
        u_grad_P_filt(:,:,:,n_data,n_filter)   = FilterFields(u_grad_P,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type);

        [du_y,du_x,du_z]                       = CentralDerivative_d1_2ndOrder(u);
        [dv_y,dv_x,dv_z]                       = CentralDerivative_d1_2ndOrder(v);
        [dw_y,dw_x,dw_z]                       = CentralDerivative_d1_2ndOrder(w);
        div_U                                  = du_x./dx + dv_y./dy + dw_z./dz;
        [sos, Beta_T, Beta_v, Beta_s, Alpha_p] = Calculate_sos(bSolver,rho,T,c_p,P,0,Fluid,Substance);
        alpha_4_filt(:,:,:,n_data,n_filter)    = FilterFields(rho.*sos.^2.*div_U,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type);

        Tau_dU         = (du_x./dx).*Tau.Tau_xx + (du_y./dy).*Tau.Tau_xy + (du_z./dz).*Tau.Tau_xz  + ...
            + (dv_x./dx).*Tau.Tau_yx + (dv_y./dy).*Tau.Tau_yy + (dv_z./dz).*Tau.Tau_yz  + ...
            + (dw_x./dx).*Tau.Tau_zx + (dw_y./dy).*Tau.Tau_zy + (dw_z./dz).*Tau.Tau_zz;
        alpha_5_unfilt = Alpha_p./(c_v.*Beta_T).*(1./rho.*(Tau_dU - q_div));
        alpha_5_filt(:,:,:,n_data,n_filter)    = FilterFields(alpha_5_unfilt,delta_filt,delta,fi(nn),x,y,z,dx,dy,dz,Filter_type);


        % Calculate bulk values, cell size and deltas
        [rho_b_filt(n_data,n_filter), mu_b_filt(n_data,n_filter), P_b_filt(n_data,n_filter), T_b_filt(n_data,n_filter), u_b_filt(n_data,n_filter), Re_b_filt(n_data,n_filter), ~, ~] = Calculate_bulk_delta_metrics(x,y,z,rho_filt(:,:,:,n_data,n_filter), mu_filt(:,:,:,n_data,n_filter), P_filt(:,:,:,n_data,n_filter), T_filt(:,:,:,n_data,n_filter), u_filt(:,:,:,n_data,n_filter), num_points_x, num_points_y, num_points_z, delta_h, dx, dy, dz);

        % First calculate average on XZ avg_rhoui, avg_rho, avg_ui
        [avg_rho_xz_filt(n_data,n_filter,:),avg_mu_xz_filt(n_data,n_filter,:),avg_u_xz_filt(n_data,n_filter,:),avg_v_xz_filt(n_data,n_filter,:),avg_w_xz_filt(n_data,n_filter,:),avg_rhou_xz_filt(n_data,n_filter,:),avg_rhov_xz_filt(n_data,n_filter,:),avg_rhow_xz_filt(n_data,n_filter,:)] = Spatial_avg_XZ(rho_filt(:,:,:,n_data,n_filter),mu_filt(:,:,:,n_data,n_filter),u_filt(:,:,:,n_data,n_filter),v_filt(:,:,:,n_data,n_filter),w_filt(:,:,:,n_data,n_filter),num_points_x,num_points_y,num_points_z);

        % Calculate favre fluctuations with spatial average
        [R_favre_uu_fluct_filt(:,:,:,n_data,n_filter),R_favre_vv_fluct_filt(:,:,:,n_data,n_filter),R_favre_ww_fluct_filt(:,:,:,n_data,n_filter)] = Calculate_favre_fluctuations_XZ(u_filt(:,:,:,n_data,n_filter),v_filt(:,:,:,n_data,n_filter),w_filt(:,:,:,n_data,n_filter),avg_rho_xz_filt(n_data,n_filter,:),avg_rhou_xz_filt(n_data,n_filter,:),avg_rhov_xz_filt(n_data,n_filter,:),avg_rhow_xz_filt(n_data,n_filter,:),num_points_y);

        % TKE
        [TKE_filt(:,:,:,n_data,n_filter), TKE_y_filt(n_data,n_filter,:),TKE_tw_filt(n_data,n_filter,:),TKE_bw_filt(n_data,n_filter,:), ...
            avg_R_favre_uu_bw_filt(n_data,n_filter,:), avg_R_favre_vv_bw_filt(n_data,n_filter,:), avg_R_favre_ww_bw_filt(n_data,n_filter,:), ...
            avg_R_favre_uu_tw_filt(n_data,n_filter,:), avg_R_favre_vv_tw_filt(n_data,n_filter,:), avg_R_favre_ww_tw_filt(n_data,n_filter,:)] = Calculate_TKE(R_favre_uu_fluct_filt(:,:,:,n_data,n_filter),R_favre_vv_fluct_filt(:,:,:,n_data,n_filter),R_favre_ww_fluct_filt(:,:,:,n_data,n_filter),num_points_x,num_points_y,num_points_z);

        % TKE total average across y
        TKE_tot_filt(nn,n_data) = sum(squeeze(TKE_y_filt(n_data,n_filter,2:end-1))'.*dy(2,2:end-1,2))/(2*delta_h);

        % Normalize wall-units
        [y_plus_bw_filt(n_data,n_filter,:),u_plus_bw_filt(n_data,n_filter,:),T_plus_bw_filt(n_data,n_filter,:),y_plus_tw_filt(n_data,n_filter,:),u_plus_tw_filt(n_data,n_filter,:),T_plus_tw_filt(n_data,n_filter,:)] = Transform_WallUnits(y,u_filt(:,:,:,n_data,n_filter),T_filt(:,:,:,n_data,n_filter),u_tau_bw,rho_bw,mu_bw,T_bw,T_tau_bw,u_tau_tw,rho_tw,mu_tw,T_tw,T_tau_tw,num_points_x,num_points_y,num_points_z,delta_h);

    

   end



end


% DNS ensembled average
if ~(n_datasets==1)
    rho_b = sum(rho_b)/n_datasets;
    mu_b  = sum(mu_b)/n_datasets;
    P_b   = sum(P_b)/n_datasets;
    T_b   = sum(T_b)/n_datasets;
    u_b   = sum(u_b)/n_datasets;
    Re_b  = sum(Re_b)/n_datasets;

    avg_rho_xz  = sum(avg_rho_xz)/n_datasets;
    avg_u_xz    = sum(avg_u_xz)/n_datasets;
    avg_v_xz    = sum(avg_v_xz)/n_datasets;
    avg_w_xz    = sum(avg_w_xz)/n_datasets;
    avg_rhou_xz = sum(avg_rhou_xz)/n_datasets;
    avg_rhov_xz = sum(avg_rhov_xz)/n_datasets;
    avg_rhow_xz = sum(avg_rhow_xz)/n_datasets;

    R_favre_uu_fluct = sum(R_favre_uu_fluct,4)/n_datasets;
    R_favre_vv_fluct = sum(R_favre_vv_fluct,4)/n_datasets;
    R_favre_ww_fluct = sum(R_favre_ww_fluct,4)/n_datasets;

    TKE    = sum(TKE,4)/n_datasets;
    TKE_y  = sum(TKE_y)/n_datasets;
    TKE_tw = sum(TKE_tw)/n_datasets;
    TKE_bw = sum(TKE_bw)/n_datasets;
    avg_R_favre_uu_bw = sum(avg_R_favre_uu_bw)/n_datasets;
    avg_R_favre_vv_bw = sum(avg_R_favre_vv_bw)/n_datasets;
    avg_R_favre_ww_bw = sum(avg_R_favre_ww_bw)/n_datasets;
    avg_R_favre_uu_tw = sum(avg_R_favre_uu_tw)/n_datasets;
    avg_R_favre_vv_tw = sum(avg_R_favre_vv_tw)/n_datasets;
    avg_R_favre_ww_tw = sum(avg_R_favre_ww_tw)/n_datasets;

    TKE_tot_filt = sum(TKE_tot_filt,2)/n_datasets;


    u_avg_dataset     = sum(u_avg_dataset,4)/n_datasets;
    v_avg_dataset     = sum(v_avg_dataset,4)/n_datasets;
    w_avg_dataset     = sum(w_avg_dataset,4)/n_datasets;
    rho_avg_dataset   = sum(rho_avg_dataset,4)/n_datasets;
    P_avg_dataset     = sum(P_avg_dataset,4)/n_datasets;
    T_avg_dataset     = sum(T_avg_dataset,4)/n_datasets;
    mu_avg_dataset    = sum(mu_avg_dataset,4)/n_datasets;
    kappa_avg_dataset = sum(kappa_avg_dataset,4)/n_datasets;
    c_p_avg_dataset   = sum(c_p_avg_dataset,4)/n_datasets;
    c_v_avg_dataset   = sum(c_v_avg_dataset,4)/n_datasets;
   

end

% Filter ensamble average
for nn = 2:length(fi)

    n_filter = nn - 1;

    u_filt_temp     = (sum(u_filt(:,:,:,:,n_filter),4))/n_datasets;
    v_filt_temp     = (sum(v_filt(:,:,:,:,n_filter),4))/n_datasets;
    w_filt_temp     = (sum(w_filt(:,:,:,:,n_filter),4))/n_datasets;
    rho_filt_temp   = (sum(rho_filt(:,:,:,:,n_filter),4))/n_datasets;
    P_filt_temp     = (sum(P_filt(:,:,:,:,n_filter),4))/n_datasets;
    T_filt_temp     = (sum(T_filt(:,:,:,:,n_filter),4))/n_datasets;
    mu_filt_temp    = (sum(mu_filt(:,:,:,:,n_filter),4))/n_datasets;
    kappa_filt_temp = (sum(kappa_filt(:,:,:,:,n_filter),4))/n_datasets;
    c_p_filt_temp   = (sum(c_p_filt(:,:,:,:,n_filter),4))/n_datasets;
    c_v_filt_temp   = (sum(c_v_filt(:,:,:,:,n_filter),4))/n_datasets;

    rho_b_filt_temp = sum(rho_b_filt(:,n_filter))/n_datasets;
    mu_b_filt_temp  = sum(mu_b_filt(:,n_filter))/n_datasets;
    P_b_filt_temp   = sum(P_b_filt(:,n_filter))/n_datasets;
    T_b_filt_temp   = sum(T_b_filt(:,n_filter))/n_datasets;
    u_b_filt_temp   = sum(u_b_filt(:,n_filter))/n_datasets;
    Re_b_filt_temp  = sum(Re_b_filt(:,n_filter))/n_datasets;

    R_favre_uu_fluct_filt_temp = sum(R_favre_uu_fluct_filt(:,:,:,:,n_filter),4)/n_datasets;
    R_favre_vv_fluct_filt_temp = sum(R_favre_vv_fluct_filt(:,:,:,:,n_filter),4)/n_datasets;
    R_favre_ww_fluct_filt_temp = sum(R_favre_ww_fluct_filt(:,:,:,:,n_filter),4)/n_datasets;

    TKE_filt_temp = sum(TKE_filt(:,:,:,:,n_filter),4)/n_datasets;

    % Other fields
    u_filt_favre_temp     = (sum(u_filt_favre(:,:,:,:,n_filter),4))/n_datasets;
    uu_filt_favre_temp    = (sum(uu_filt_favre(:,:,:,:,n_filter),4))/n_datasets;
    v_filt_favre_temp     = (sum(v_filt_favre(:,:,:,:,n_filter),4))/n_datasets;
    vv_filt_favre_temp    = (sum(vv_filt_favre(:,:,:,:,n_filter),4))/n_datasets;
    w_filt_favre_temp     = (sum(w_filt_favre(:,:,:,:,n_filter),4))/n_datasets;
    ww_filt_favre_temp    = (sum(ww_filt_favre(:,:,:,:,n_filter),4))/n_datasets;
    uv_filt_favre_temp    = (sum(uv_filt_favre(:,:,:,:,n_filter),4))/n_datasets;
    uw_filt_favre_temp    = (sum(uw_filt_favre(:,:,:,:,n_filter),4))/n_datasets;
    vw_filt_favre_temp    = (sum(vw_filt_favre(:,:,:,:,n_filter),4))/n_datasets;
    rhou_visc_filt_temp   = (sum(rhou_visc_filt(:,:,:,:,n_filter),4))/n_datasets;
    rhov_visc_filt_temp   = (sum(rhov_visc_filt(:,:,:,:,n_filter),4))/n_datasets;
    rhow_visc_filt_temp   = (sum(rhow_visc_filt(:,:,:,:,n_filter),4))/n_datasets;
    u_grad_P_filt_temp    = (sum(u_grad_P_filt(:,:,:,:,n_filter),4))/n_datasets;
    alpha_4_filt_temp     = (sum(alpha_4_filt(:,:,:,:,n_filter),4))/n_datasets;
    alpha_5_filt_temp     = (sum(alpha_5_filt(:,:,:,:,n_filter),4))/n_datasets;



    if ~(n_datasets==1)
        avg_rho_xz_filt_temp   = squeeze(sum(avg_rho_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_mu_xz_filt_temp    = squeeze(sum(avg_mu_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_u_xz_filt_temp     = squeeze(sum(avg_u_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_v_xz_filt_temp     = squeeze(sum(avg_v_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_w_xz_filt_temp     = squeeze(sum(avg_w_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_rhou_xz_filt_temp  = squeeze(sum(avg_rhou_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_rhov_xz_filt_temp  = squeeze(sum(avg_rhov_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_rhow_xz_filt_temp  = squeeze(sum(avg_rhow_xz_filt(:,n_filter,:)))'/n_datasets;

        TKE_y_filt_temp  = squeeze(sum(TKE_y_filt(:,n_filter,:)))'/n_datasets;
        TKE_tw_filt_temp = squeeze(sum(TKE_tw_filt(:,n_filter,:)))'/n_datasets;
        TKE_bw_filt_temp = squeeze(sum(TKE_bw_filt(:,n_filter,:)))'/n_datasets;

        avg_R_favre_uu_bw_filt_temp = squeeze(sum(avg_R_favre_uu_bw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_vv_bw_filt_temp = squeeze(sum(avg_R_favre_vv_bw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_ww_bw_filt_temp = squeeze(sum(avg_R_favre_ww_bw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_uu_tw_filt_temp = squeeze(sum(avg_R_favre_uu_tw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_vv_tw_filt_temp = squeeze(sum(avg_R_favre_vv_tw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_ww_tw_filt_temp = squeeze(sum(avg_R_favre_ww_tw_filt(:,n_filter,:)))'/n_datasets;


        y_plus_bw_filt_temp = squeeze(sum(y_plus_bw_filt(:,n_filter,:),1))'/n_datasets;
        u_plus_bw_filt_temp = squeeze(sum(u_plus_bw_filt(:,n_filter,:),1))'/n_datasets;
        T_plus_bw_filt_temp = squeeze(sum(T_plus_bw_filt(:,n_filter,:),1))'/n_datasets;
        y_plus_tw_filt_temp = squeeze(sum(y_plus_tw_filt(:,n_filter,:),1))'/n_datasets;
        u_plus_tw_filt_temp = squeeze(sum(u_plus_tw_filt(:,n_filter,:),1))'/n_datasets;
        T_plus_tw_filt_temp = squeeze(sum(T_plus_tw_filt(:,n_filter,:),1))'/n_datasets;

    else

        avg_rho_xz_filt_temp   = squeeze((avg_rho_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_mu_xz_filt_temp    = squeeze((avg_mu_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_u_xz_filt_temp     = squeeze((avg_u_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_v_xz_filt_temp     = squeeze((avg_v_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_w_xz_filt_temp     = squeeze((avg_w_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_rhou_xz_filt_temp  = squeeze((avg_rhou_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_rhov_xz_filt_temp  = squeeze((avg_rhov_xz_filt(:,n_filter,:)))'/n_datasets;
        avg_rhow_xz_filt_temp  = squeeze((avg_rhow_xz_filt(:,n_filter,:)))'/n_datasets;

        TKE_y_filt_temp  = squeeze((TKE_y_filt(:,n_filter,:)))'/n_datasets;
        TKE_tw_filt_temp = squeeze((TKE_tw_filt(:,n_filter,:)))'/n_datasets;
        TKE_bw_filt_temp = squeeze((TKE_bw_filt(:,n_filter,:)))'/n_datasets;

        avg_R_favre_uu_bw_filt_temp = squeeze((avg_R_favre_uu_bw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_vv_bw_filt_temp = squeeze((avg_R_favre_vv_bw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_ww_bw_filt_temp = squeeze((avg_R_favre_ww_bw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_uu_tw_filt_temp = squeeze((avg_R_favre_uu_tw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_vv_tw_filt_temp = squeeze((avg_R_favre_vv_tw_filt(:,n_filter,:)))'/n_datasets;
        avg_R_favre_ww_tw_filt_temp = squeeze((avg_R_favre_ww_tw_filt(:,n_filter,:)))'/n_datasets;


        y_plus_bw_filt_temp = squeeze((y_plus_bw_filt(:,n_filter,:)))'/n_datasets;
        u_plus_bw_filt_temp = squeeze((u_plus_bw_filt(:,n_filter,:)))'/n_datasets;
        T_plus_bw_filt_temp = squeeze((T_plus_bw_filt(:,n_filter,:)))'/n_datasets;
        y_plus_tw_filt_temp = squeeze((y_plus_tw_filt(:,n_filter,:)))'/n_datasets;
        u_plus_tw_filt_temp = squeeze((u_plus_tw_filt(:,n_filter,:)))'/n_datasets;
        T_plus_tw_filt_temp = squeeze((T_plus_tw_filt(:,n_filter,:)))'/n_datasets;



    end


    % Store values
    assignin('base',strcat('TKE_filt',num2str(fi(nn)),'xDelta'), TKE_filt_temp)
    assignin('base',strcat('TKE_y_filt',num2str(fi(nn)),'xDelta'), TKE_y_filt_temp)
    assignin('base',strcat('TKE_tw_filt',num2str(fi(nn)),'xDelta'), TKE_tw_filt_temp)
    assignin('base',strcat('TKE_bw_filt',num2str(fi(nn)),'xDelta'), TKE_bw_filt_temp)

    assignin('base',strcat('y_plus_bw_filt',num2str(fi(nn)),'xDelta'), y_plus_bw_filt_temp)
    assignin('base',strcat('u_plus_bw_filt',num2str(fi(nn)),'xDelta'), u_plus_bw_filt_temp)
    assignin('base',strcat('T_plus_bw_filt',num2str(fi(nn)),'xDelta'), T_plus_bw_filt_temp)
    assignin('base',strcat('y_plus_tw_filt',num2str(fi(nn)),'xDelta'), y_plus_tw_filt_temp)
    assignin('base',strcat('u_plus_tw_filt',num2str(fi(nn)),'xDelta'), u_plus_tw_filt_temp)
    assignin('base',strcat('T_plus_tw_filt',num2str(fi(nn)),'xDelta'), T_plus_tw_filt_temp)
    assignin('base',strcat('avg_R_favre_uu_bw_filt',num2str(fi(nn)),'xDelta'), avg_R_favre_uu_bw_filt_temp)
    assignin('base',strcat('avg_R_favre_vv_bw_filt',num2str(fi(nn)),'xDelta'), avg_R_favre_vv_bw_filt_temp)
    assignin('base',strcat('avg_R_favre_ww_bw_filt',num2str(fi(nn)),'xDelta'), avg_R_favre_ww_bw_filt_temp)
    assignin('base',strcat('avg_R_favre_uu_tw_filt',num2str(fi(nn)),'xDelta'), avg_R_favre_uu_tw_filt_temp)
    assignin('base',strcat('avg_R_favre_vv_tw_filt',num2str(fi(nn)),'xDelta'), avg_R_favre_vv_tw_filt_temp)
    assignin('base',strcat('avg_R_favre_ww_tw_filt',num2str(fi(nn)),'xDelta'), avg_R_favre_ww_tw_filt_temp)

    assignin('base',strcat('avg_rho_xz_filt',num2str(fi(nn)),'xDelta'), avg_rho_xz_filt_temp)
    assignin('base',strcat('avg_mu_xz_filt',num2str(fi(nn)),'xDelta'), avg_mu_xz_filt_temp)
    assignin('base',strcat('avg_u_xz_filt',num2str(fi(nn)),'xDelta'), avg_u_xz_filt_temp)
    assignin('base',strcat('avg_v_xz_filt',num2str(fi(nn)),'xDelta'), avg_v_xz_filt_temp)
    assignin('base',strcat('avg_w_xz_filt',num2str(fi(nn)),'xDelta'), avg_w_xz_filt_temp)
    assignin('base',strcat('avg_rhou_xz_filt',num2str(fi(nn)),'xDelta'), avg_rhou_xz_filt_temp)
    assignin('base',strcat('avg_rhov_xz_filt',num2str(fi(nn)),'xDelta'), avg_rhov_xz_filt_temp)
    assignin('base',strcat('avg_rhow_xz_filt',num2str(fi(nn)),'xDelta'), avg_rhow_xz_filt_temp)

    assignin('base',strcat('u_b_filt',num2str(fi(nn)),'xDelta'), u_b_filt_temp)

    % 3D filtered data
    assignin('base',strcat('u_filt',num2str(fi(nn)),'xDelta'), u_filt_temp)
    assignin('base',strcat('v_filt',num2str(fi(nn)),'xDelta'), v_filt_temp)
    assignin('base',strcat('w_filt',num2str(fi(nn)),'xDelta'), w_filt_temp)
    assignin('base',strcat('rho_filt',num2str(fi(nn)),'xDelta'), rho_filt_temp)
    assignin('base',strcat('P_filt',num2str(fi(nn)),'xDelta'), P_filt_temp)
    assignin('base',strcat('T_filt',num2str(fi(nn)),'xDelta'), T_filt_temp)
    assignin('base',strcat('mu_filt',num2str(fi(nn)),'xDelta'), mu_filt_temp)
    assignin('base',strcat('kappa_filt',num2str(fi(nn)),'xDelta'), kappa_filt_temp)
    assignin('base',strcat('c_p_filt',num2str(fi(nn)),'xDelta'), c_p_filt_temp)
    assignin('base',strcat('c_v_filt',num2str(fi(nn)),'xDelta'), c_v_filt_temp)

    assignin('base',strcat('R_favre_uu_fluct_filt',num2str(fi(nn)),'xDelta'), R_favre_uu_fluct_filt_temp)
    assignin('base',strcat('R_favre_vv_fluct_filt',num2str(fi(nn)),'xDelta'), R_favre_vv_fluct_filt_temp)
    assignin('base',strcat('R_favre_ww_fluct_filt',num2str(fi(nn)),'xDelta'), R_favre_ww_fluct_filt_temp)


    % Other snapshot average filtered fields
    assignin('base',strcat('u_filt_favre',num2str(fi(nn)),'xDelta'), u_filt_favre_temp)
    assignin('base',strcat('uu_filt_favre',num2str(fi(nn)),'xDelta'), uu_filt_favre_temp)
    assignin('base',strcat('v_filt_favre',num2str(fi(nn)),'xDelta'), v_filt_favre_temp)
    assignin('base',strcat('vv_filt_favre',num2str(fi(nn)),'xDelta'), vv_filt_favre_temp)
    assignin('base',strcat('w_filt_favre',num2str(fi(nn)),'xDelta'), w_filt_favre_temp)
    assignin('base',strcat('ww_filt_favre',num2str(fi(nn)),'xDelta'), ww_filt_favre_temp)
    assignin('base',strcat('uv_filt_favre',num2str(fi(nn)),'xDelta'), uv_filt_favre_temp)
    assignin('base',strcat('uw_filt_favre',num2str(fi(nn)),'xDelta'), uw_filt_favre_temp)
    assignin('base',strcat('vw_filt_favre',num2str(fi(nn)),'xDelta'), vw_filt_favre_temp)
    assignin('base',strcat('rhou_visc_filt',num2str(fi(nn)),'xDelta'), rhou_visc_filt_temp)
    assignin('base',strcat('rhov_visc_filt',num2str(fi(nn)),'xDelta'), rhov_visc_filt_temp)
    assignin('base',strcat('rhow_visc_filt',num2str(fi(nn)),'xDelta'), rhow_visc_filt_temp)
    assignin('base',strcat('u_grad_P_filt',num2str(fi(nn)),'xDelta'), u_grad_P_filt_temp)
    assignin('base',strcat('alpha_4_filt',num2str(fi(nn)),'xDelta'), alpha_4_filt_temp)
    assignin('base',strcat('alpha_5_filt',num2str(fi(nn)),'xDelta'), alpha_5_filt_temp)

end



%% Save DNS and filtered results
if bSave

    %% Save ensemble-average (XZ) into .csv file
    name_file_out = strcat('HP_','fi_2_10_single', Filter_type);
    DataOutput_LES(y(2,:,2), avg_u_xz, u_b, y_plus_bw,u_plus_bw, T_plus_bw, y_plus_tw,u_plus_tw, T_plus_tw, avg_R_favre_uu_bw, avg_R_favre_vv_bw, avg_R_favre_ww_bw, avg_R_favre_uu_tw, avg_R_favre_vv_tw, avg_R_favre_ww_tw, TKE_y, TKE_bw, TKE_tw, name_file_out, ...
        avg_u_xz_filt2xDelta, u_b_filt2xDelta, y_plus_bw_filt2xDelta, u_plus_bw_filt2xDelta, T_plus_bw_filt2xDelta, y_plus_tw_filt2xDelta, u_plus_tw_filt2xDelta, T_plus_tw_filt2xDelta, avg_R_favre_uu_bw_filt2xDelta, avg_R_favre_vv_bw_filt2xDelta, avg_R_favre_ww_bw_filt2xDelta, avg_R_favre_uu_tw_filt2xDelta, avg_R_favre_vv_tw_filt2xDelta, avg_R_favre_ww_tw_filt2xDelta, TKE_y_filt2xDelta, TKE_bw_filt2xDelta, TKE_tw_filt2xDelta,...
        avg_u_xz_filt4xDelta, u_b_filt4xDelta, y_plus_bw_filt4xDelta, u_plus_bw_filt4xDelta, T_plus_bw_filt4xDelta, y_plus_tw_filt4xDelta, u_plus_tw_filt4xDelta, T_plus_tw_filt4xDelta, avg_R_favre_uu_bw_filt4xDelta, avg_R_favre_vv_bw_filt4xDelta, avg_R_favre_ww_bw_filt4xDelta, avg_R_favre_uu_tw_filt4xDelta, avg_R_favre_vv_tw_filt4xDelta, avg_R_favre_ww_tw_filt4xDelta, TKE_y_filt4xDelta, TKE_bw_filt4xDelta, TKE_tw_filt4xDelta,...
        avg_u_xz_filt6xDelta, u_b_filt6xDelta, y_plus_bw_filt6xDelta, u_plus_bw_filt6xDelta, T_plus_bw_filt6xDelta, y_plus_tw_filt6xDelta, u_plus_tw_filt6xDelta, T_plus_tw_filt6xDelta, avg_R_favre_uu_bw_filt6xDelta, avg_R_favre_vv_bw_filt6xDelta, avg_R_favre_ww_bw_filt6xDelta, avg_R_favre_uu_tw_filt6xDelta, avg_R_favre_vv_tw_filt6xDelta, avg_R_favre_ww_tw_filt6xDelta, TKE_y_filt6xDelta, TKE_bw_filt6xDelta, TKE_tw_filt6xDelta,...
        avg_u_xz_filt8xDelta, u_b_filt8xDelta, y_plus_bw_filt8xDelta, u_plus_bw_filt8xDelta, T_plus_bw_filt8xDelta, y_plus_tw_filt8xDelta, u_plus_tw_filt8xDelta, T_plus_tw_filt8xDelta, avg_R_favre_uu_bw_filt8xDelta, avg_R_favre_vv_bw_filt8xDelta, avg_R_favre_ww_bw_filt8xDelta, avg_R_favre_uu_tw_filt8xDelta, avg_R_favre_vv_tw_filt8xDelta, avg_R_favre_ww_tw_filt8xDelta, TKE_y_filt8xDelta, TKE_bw_filt8xDelta, TKE_tw_filt8xDelta,...
        avg_u_xz_filt10xDelta, u_b_filt10xDelta, y_plus_bw_filt10xDelta, u_plus_bw_filt10xDelta, T_plus_bw_filt10xDelta, y_plus_tw_filt10xDelta, u_plus_tw_filt10xDelta, T_plus_tw_filt10xDelta, avg_R_favre_uu_bw_filt10xDelta, avg_R_favre_vv_bw_filt10xDelta, avg_R_favre_ww_bw_filt10xDelta, avg_R_favre_uu_tw_filt10xDelta, avg_R_favre_vv_tw_filt10xDelta, avg_R_favre_ww_tw_filt10xDelta, TKE_y_filt10xDelta, TKE_bw_filt10xDelta, TKE_tw_filt10xDelta);

    %% Export all filtered fields in .mat file for postprocess
    name_file_out_filt_dataset = strcat("Filtered_DNS_ensemble_average_fi_single");
    DataOutput_EnsembleAverage_Filt;

end

if bPlot

    %% Plots based on latest snapshot
    % Filtered wall units based on DNS tau values
    Plot_filter;

end

toc
