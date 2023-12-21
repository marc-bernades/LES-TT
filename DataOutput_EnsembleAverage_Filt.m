%% Export ensemble-average filtered DNS at 2 and 4 Delta

% Save to mat file
Headings_Main        = {'x', 'y', 'z', 'delta', 'dx', 'dy', 'dz'}; 
Heading_Average_Db   = {'u_avg_dataset','v_avg_dataset','w_avg_dataset','rho_avg_dataset','P_avg_dataset','T_avg_dataset','mu_avg_dataset','kappa_avg_dataset','c_p_avg_dataset','c_v_avg_dataset'};
Main                 = {x, y, z, delta, dx, dy, dz};
Main_Average_Dataset = {u_avg_dataset,v_avg_dataset,w_avg_dataset,rho_avg_dataset,P_avg_dataset,T_avg_dataset,mu_avg_dataset,kappa_avg_dataset,c_p_avg_dataset,c_v_avg_dataset};

% Ensemble headings and data
Data_output   = [Headings_Main Heading_Average_Db; [Main Main_Average_Dataset]];


Headings_Filt = cell(length(fi) - 1);
Data_Filt     = cell(length(fi) - 1);

for nn = 2:length(fi)


Headings_Filt{nn-1}     = {"rho_filt" + num2str(fi(nn)) + "xDelta", "u_filt" + num2str(fi(nn)) + "xDelta", "v_filt" + num2str(fi(nn)) + "xDelta", "w_filt" + num2str(fi(nn)) + "xDelta", "P_filt" + num2str(fi(nn)) + "xDelta", "T_filt" + num2str(fi(nn)) + "xDelta", "mu_filt" + num2str(fi(nn)) + "xDelta", "kappa_filt" + num2str(fi(nn)) + "xDelta", "c_p_filt" + num2str(fi(nn)) + "xDelta", "c_v_filt" + num2str(fi(nn)) + "xDelta",...
    "u_filt_favre"  + num2str(fi(nn)) + "xDelta", "uu_filt_favre" + num2str(fi(nn)) + "xDelta", "v_filt_favre" + num2str(fi(nn)) + "xDelta", "vv_filt_favre" + num2str(fi(nn)) + "xDelta","w_filt_favre" + num2str(fi(nn)) + "xDelta","ww_filt_favre" + num2str(fi(nn)) + "xDelta"...
    "uv_filt_favre" + num2str(fi(nn)) + "xDelta", "uw_filt_favre" + num2str(fi(nn)) + "xDelta", "vw_filt_favre" + num2str(fi(nn)) + "xDelta", "rhou_visc_filt" + num2str(fi(nn)) + "xDelta", "rhov_visc_filt" + num2str(fi(nn)) + "xDelta","rhow_visc_filt" + num2str(fi(nn)) + "xDelta", ...
    "u_grad_P_filt" + num2str(fi(nn)) + "xDelta", "alpha_4_filt" + num2str(fi(nn)) + "xDelta", "alpha_5_filt" + num2str(fi(nn)) + "xDelta"};


Data_Filt{nn-1}     = {eval("rho_filt" + num2str(fi(nn)) + "xDelta"), eval("u_filt" + num2str(fi(nn)) + "xDelta"), eval("v_filt" + num2str(fi(nn)) + "xDelta"), eval("w_filt" + num2str(fi(nn)) + "xDelta"), eval("P_filt" + num2str(fi(nn)) + "xDelta"), eval("T_filt" + num2str(fi(nn)) + "xDelta"), eval("mu_filt" + num2str(fi(nn)) + "xDelta"), eval("kappa_filt" + num2str(fi(nn)) + "xDelta"), eval("c_p_filt" + num2str(fi(nn)) + "xDelta"), eval("c_v_filt" + num2str(fi(nn)) + "xDelta"),...
eval("u_filt_favre"  + num2str(fi(nn)) + "xDelta"), eval("uu_filt_favre" + num2str(fi(nn)) + "xDelta"), eval("v_filt_favre" + num2str(fi(nn)) + "xDelta"), eval("vv_filt_favre" + num2str(fi(nn)) + "xDelta"), eval("w_filt_favre" + num2str(fi(nn)) + "xDelta"), eval("ww_filt_favre" + num2str(fi(nn)) + "xDelta"),...
eval("uv_filt_favre" + num2str(fi(nn)) + "xDelta"), eval("uw_filt_favre" + num2str(fi(nn)) + "xDelta"), eval("vw_filt_favre" + num2str(fi(nn)) + "xDelta"), eval("rhou_visc_filt" + num2str(fi(nn)) + "xDelta"), eval("rhov_visc_filt" + num2str(fi(nn)) + "xDelta"), eval("rhow_visc_filt" + num2str(fi(nn)) + "xDelta"), ...
eval("u_grad_P_filt" + num2str(fi(nn)) + "xDelta"), eval("alpha_4_filt" + num2str(fi(nn)) + "xDelta"), eval("alpha_5_filt" + num2str(fi(nn)) + "xDelta")};

Data_filt_append = [Headings_Filt{nn-1}; Data_Filt{nn-1}];
Data_output      = [Data_output Data_filt_append];

end

% Convert cell to a table and use first row as variable names
save(["Data/" name_file_out_filt_dataset ".mat"], 'Data_output', '-v7.3');

