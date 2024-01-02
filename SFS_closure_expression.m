%% Closure expressions

%% Alpha 2
% Reference
rhou_visc_filt_avg       = Spatial_avg_XZ_var(rhou_visc_filt,num_points_x,num_points_y,num_points_z, fi, idx);
rhov_visc_filt_avg       = Spatial_avg_XZ_var(rhov_visc_filt,num_points_x,num_points_y,num_points_z, fi, idx);
rhow_visc_filt_avg       = Spatial_avg_XZ_var(rhow_visc_filt,num_points_x,num_points_y,num_points_z, fi, idx);

%% Alpha 2
% Method 1 - Gradient
[drhou_visc_LES_y,drhou_visc_LES_x,drhou_visc_LES_z] = CentralDerivative_d1_2ndOrder(rhou_visc_LES);
[drhov_visc_LES_y,drhov_visc_LES_x,drhov_visc_LES_z] = CentralDerivative_d1_2ndOrder(rhov_visc_LES);
[drhow_visc_LES_y,drhow_visc_LES_x,drhow_visc_LES_z] = CentralDerivative_d1_2ndOrder(rhow_visc_LES);

alpha_2_closure_grad{1}       = (drhou_visc_LES_x./dx + drhou_visc_LES_y./dy + drhou_visc_LES_z./dz);
alpha_2_closure_grad{2}       = (drhov_visc_LES_x./dx + drhov_visc_LES_y./dy + drhov_visc_LES_z./dz);
alpha_2_closure_grad{3}       = (drhow_visc_LES_x./dx + drhow_visc_LES_y./dy + drhow_visc_LES_z./dz);

v_SFS_2_sim = 1E-5;

alpha_2_closure_grad_avg{1}   = v_SFS_2_sim*Spatial_avg_XZ_var(alpha_2_closure_grad{1},num_points_x,num_points_y,num_points_z, fi, idx);
alpha_2_closure_grad_avg{2}   = v_SFS_2_sim*Spatial_avg_XZ_var(alpha_2_closure_grad{2},num_points_x,num_points_y,num_points_z, fi, idx);
alpha_2_closure_grad_avg{3}   = v_SFS_2_sim*Spatial_avg_XZ_var(alpha_2_closure_grad{3},num_points_x,num_points_y,num_points_z, fi, idx);

% % Momentum X
% figure; hold on; box on
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhou_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhou_visc_filt_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_grad_avg{1}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.8500 0.3250 0.0980])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_grad_avg{1}(idx:end-idx-1) + rhou_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-.','color','k')
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_avg{1}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])
% 
% legend('LES','filt','Closure','LES + Closure','alpha2')
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel(strcat('${\langle || (\overline{\rho} \widetilde{u} / t)^\ast || \rangle}$', ' Decomposition'),'interpreter','latex')
% pbaspect([1.2 1 1])
% % legend('Location','northwest','box','off','NumColumns', 2)
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',14)
% % xlim([0 2]); ylim([-0.1 0.1])
% symlog(gca,'y',-4); grid off;
% 
% % Momentum Y
% figure; hold on; box on
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhov_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhov_visc_filt_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_grad_avg{2}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.8500 0.3250 0.0980])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_grad_avg{2}(idx:end-idx-1) + rhov_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-.','color','k')
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_avg{2}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])
% 
% legend('LES','filt','Closure','LES + Closure','alpha2')
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel(strcat('${\langle || (\overline{\rho} \widetilde{u} / t)^\ast || \rangle}$', ' Decomposition'),'interpreter','latex')
% pbaspect([1.2 1 1])
% % legend('Location','northwest','box','off','NumColumns', 2)
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',14)
% % xlim([0 2]); ylim([-0.1 0.1])
% symlog(gca,'y',-4); grid off;
% 
% % Momentum Z
% figure; hold on; box on
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhow_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhow_visc_filt_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_grad_avg{3}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.8500 0.3250 0.0980])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_grad_avg{3}(idx:end-idx-1) + rhow_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-.','color','k')
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_avg{3}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])
% 
% legend('LES','filt','Closure','LES + Closure','alpha2')
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel(strcat('${\langle || (\overline{\rho} \widetilde{u} / t)^\ast || \rangle}$', ' Decomposition'),'interpreter','latex')
% pbaspect([1.2 1 1])
% % legend('Location','northwest','box','off','NumColumns', 2)
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',14)
% % xlim([0 2]); ylim([-0.1 0.1])
% symlog(gca,'y',-4); grid off;


% Method 2 - Scale similarity
u_filt_favre_filt    = FilterFields(rho_filt.*u_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_filt;
v_filt_favre_filt    = FilterFields(rho_filt.*v_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_filt;
w_filt_favre_filt    = FilterFields(rho_filt.*w_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_filt;
mu_filt_filt         = FilterFields(mu_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
kappa_filt_filt      = FilterFields(kappa_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
T_filt_filt          = FilterFields(T_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);


[rhou_visc_LES_filt, rhov_visc_LES_filt, rhow_visc_LES_filt, rhoE_visc_LES_filt, Tau_LES_filt, q_div_LES_filt] = Calculate_Viscous_LES(u_filt_favre_filt,v_filt_favre_filt,w_filt_favre_filt,mu_filt_filt,kappa_filt_filt,T_filt_filt,dx,dy,dz,1,1);
rhou_visc_filt_filt        = FilterFields(rhou_visc_LES,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rhov_visc_filt_filt        = FilterFields(rhov_visc_LES,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
rhow_visc_filt_filt        = FilterFields(rhow_visc_LES,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);


alpha_2_closure_sim{1}     = rhou_visc_filt_filt - rhou_visc_LES_filt;
alpha_2_closure_sim{2}     = rhov_visc_filt_filt - rhov_visc_LES_filt;
alpha_2_closure_sim{3}     = rhow_visc_filt_filt - rhow_visc_LES_filt;

v_SFS_2_sim = 40;
alpha_2_closure_sim_u_avg = Spatial_avg_XZ_var(alpha_2_closure_sim{1},num_points_x,num_points_y,num_points_z, fi, idx);
v_SFS_2_sim    = 2;
% v_SFS_2_sim    = 1/mean((rhou_visc_filt_avg(idx:end-idx-1) - rhou_visc_LES_avg(idx:end-idx-1))./alpha_2_closure_sim_u_avg(idx:end-idx-1));


alpha_2_closure_sim_avg{1}  = v_SFS_2_sim*Spatial_avg_XZ_var(alpha_2_closure_sim{1},num_points_x,num_points_y,num_points_z, fi, idx);
alpha_2_closure_sim_avg{2}  = v_SFS_2_sim*Spatial_avg_XZ_var(alpha_2_closure_sim{2},num_points_x,num_points_y,num_points_z, fi, idx);
alpha_2_closure_sim_avg{3}  = v_SFS_2_sim*Spatial_avg_XZ_var(alpha_2_closure_sim{3},num_points_x,num_points_y,num_points_z, fi, idx);


% Momentum X

figure; hold on; box on
plot(y(2,idx:end-idx-1,2)/delta_h,(rhou_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
plot(y(2,idx:end-idx-1,2)/delta_h,(rhou_visc_filt_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_sim_avg{1}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color','k')
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_sim_avg{1}(idx:end-idx-1) + rhou_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-.','color',[0.4660 0.6740 0.1880])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_avg{1}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])

legend('LES','filt','Closure','LES + Closure')
xlabel('${y/\delta}$','interpreter','latex')
ylabel(strcat('${\langle (\partial (\overline{\rho} \widetilde{u}) / \partial t)^\ast \rangle}$', ' Decomposition'),'interpreter','latex')
legend([{'$(\nabla \cdot \mathbf{\breve{\sigma}})^\ast$'},{'$(\overline{\nabla \cdot \mathbf{{\sigma}}})^\ast$'}, ...
    {'${{\alpha_2}^{SFS}}^\ast$'},{'$(\nabla \cdot \mathbf{\breve{\sigma}} + {{\alpha_2}^{SFS}})^\ast$'}],'interpreter','latex','location','best')
% legend([{'$(\partial \overline{\rho} \widetilde{u} \widetilde{u} / \partial x)^\ast$'},{'$(\partial \breve{\sigma_x} / \partial x)^\ast$'},{'$(\partial \overline{P} / \partial x)^\ast$'},{'$\alpha_1^\ast$'},{'$\alpha_2^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.2 1 1])
legend('Location','northwest','box','off','NumColumns', 2)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]); ylim([-0.01 0.01])
symlog(gca,'y',-4); grid off;
exportgraphics(gca,strcat('Figures/SFS_Alpha2_', num2str(fi), 'xDelta', '.jpeg'),'Resolution',300)
exportgraphics(gca,strcat('Figures/SFS_Alpha2_', num2str(fi), 'xDelta', '.png'),'Resolution',300)


% % Momentum Y
% figure; hold on; box on
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhov_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhov_visc_filt_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_sim_avg{2}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.8500 0.3250 0.0980])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_sim_avg{2}(idx:end-idx-1) + rhov_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-.','color','k')
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_avg{2}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])
% 
% legend('LES','filt','Closure','LES + Closure','alpha2')
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel(strcat('${\langle || (\overline{\rho} \widetilde{u} / t)^\ast || \rangle}$', ' Decomposition'),'interpreter','latex')
% pbaspect([1.2 1 1])
% % legend('Location','northwest','box','off','NumColumns', 2)
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',14)
% xlim([0 2]); ylim([-0.0001 0.0001])
% symlog(gca,'y',-6); grid off;

% % Momentum Z
% figure; hold on; box on
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhow_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
% plot(y(2,idx:end-idx-1,2)/delta_h,(rhow_visc_filt_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_sim_avg{3}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.8500 0.3250 0.0980])
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_closure_sim_avg{3}(idx:end-idx-1) + rhow_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-.','color','k')
% plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_avg{3}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])
% 
% legend('LES','filt','Closure','LES + Closure','alpha2')
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel(strcat('${\langle || (\overline{\rho} \widetilde{u} / t)^\ast || \rangle}$', ' Decomposition'),'interpreter','latex')
% pbaspect([1.2 1 1])
% % legend('Location','northwest','box','off','NumColumns', 2)
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',14)
% xlim([0 2]); ylim([-0.001 0.001])
% symlog(gca,'y',-5); grid off;




%% Alpha 4
alpha_4_SFS_filt      = FilterFields(rho_filt.*sos_filt.^2.*divU_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);

% Method 1 - Scale similarity
sos_filt_filt         = FilterFields(sos_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
% u_filt_favre_favre    = FilterFields(rho_filt.*u_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_filt;
% v_filt_favre_favre    = FilterFields(rho_filt.*v_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_filt;
% w_filt_favre_favre    = FilterFields(rho_filt.*w_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_filt;
u_filt_favre_favre    = FilterFields(u_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
v_filt_favre_favre    = FilterFields(v_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
w_filt_favre_favre    = FilterFields(w_filt_favre,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);

[du_y_filt_favre_favre,du_x_filt_favre_favre,du_z_filt_favre_favre] = CentralDerivative_d1_2ndOrder(u_filt_favre_favre);
[dv_y_filt_favre_favre,dv_x_filt_favre_favre,dv_z_filt_favre_favre] = CentralDerivative_d1_2ndOrder(v_filt_favre_favre);
[dw_y_filt_favre_favre,dw_x_filt_favre_favre,dw_z_filt_favre_favre] = CentralDerivative_d1_2ndOrder(w_filt_favre_favre);

% Strain tensor
divU_filt_favre_favre = du_x_filt_favre_favre./dx + dv_y_filt_favre_favre./dy + dw_z_filt_favre_favre./dz;
alpha_4_closure       = rho_filt_filt.*sos_filt_filt.^2.*divU_filt_favre_favre - alpha_4_SFS_filt;
alpha_4_closure_avg   = Spatial_avg_XZ_var(alpha_4_closure,num_points_x,num_points_y,num_points_z, fi, idx);
alpha_4_filt_avg      = Spatial_avg_XZ_var(alpha_4_filt,num_points_x,num_points_y,num_points_z, fi, idx);

% Method 1 - Germano identity version
alpha_4_SFS_FG     = FilterFields(rho_filt.*sos_filt.^2.*divU_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);

sos_filt_FG        = FilterFields(sos_filt,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
% u_filt_favre_FG    = FilterFields(rho_filt.*u_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_FG_filt;
% v_filt_favre_FG    = FilterFields(rho_filt.*v_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_FG_filt;
% w_filt_favre_FG    = FilterFields(rho_filt.*w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type)./rho_filt_FG_filt;
u_filt_favre_FG    = FilterFields(u_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
v_filt_favre_FG    = FilterFields(v_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);
w_filt_favre_FG    = FilterFields(w_filt_favre,2*delta_filt,delta,2*fi,x,y,z,dx,dy,dz,Filter_type);

[du_y_filt_favre_FG,du_x_filt_favre_FG,du_z_filt_favre_FG] = CentralDerivative_d1_2ndOrder(u_filt_favre_FG);
[dv_y_filt_favre_FG,dv_x_filt_favre_FG,dv_z_filt_favre_FG] = CentralDerivative_d1_2ndOrder(v_filt_favre_FG);
[dw_y_filt_favre_FG,dw_x_filt_favre_FG,dw_z_filt_favre_FG] = CentralDerivative_d1_2ndOrder(w_filt_favre_FG);

% Strain tensor
divU_filt_favre_FG       = du_x_filt_favre_FG./dx + dv_y_filt_favre_FG./dy + dw_z_filt_favre_FG./dz;
alpha_4_closure_FG       = rho_filt_FG.*sos_filt_FG.^2.*divU_filt_favre_FG - alpha_4_SFS_FG;
alpha_4_closure_FG_avg   = Spatial_avg_XZ_var(alpha_4_closure_FG,num_points_x,num_points_y,num_points_z, fi, idx);

% Coefficients
v_SFS    = -1/mean((alpha_4_filt_avg(idx:end-idx-1) - alpha_4_LES_avg(idx:end-idx-1))./alpha_4_closure_avg(idx:end-idx-1));
v_SFS_FG = -1/mean((alpha_4_filt_avg(idx:end-idx-1) - alpha_4_LES_avg(idx:end-idx-1))./alpha_4_closure_FG_avg(idx:end-idx-1));
v_SFS = -10;

% Figure Pressure Method 1
figure; hold on; box on
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_4_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_4_filt_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_4_closure_FG_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','--','color','k')
plot(y(2,idx:end-idx-1,2)/delta_h,(-1/(v_SFS)*alpha_4_closure_avg(idx:end-idx-1) + alpha_4_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-.','color',[0.4660 0.6740 0.1880])

xlabel('${y/\delta}$','interpreter','latex')
ylabel(strcat('${\langle (\partial P / \partial t)^\ast \rangle}$', 'Decomposition'),'interpreter','latex')
legend([{'$(\overline{\rho} \breve{c}^2 \nabla \cdot \widetilde{u})^\ast$'},{'$(\overline{{\rho} {c}^2 \nabla \cdot {u}})^\ast$'},...
    {'$({\alpha_4}^{SFS})^\ast$'},{'$(\overline{\rho} \breve{c}^2 \nabla \cdot \widetilde{u} + {\alpha_4}^{SFS})^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.2 1 1])
legend('Location','northeast','box','off','NumColumns', 2)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]); ylim([-1000 1000])
% ylim([10^-2 10^4])
% set(gca,'YTick',10.^(-2:2:4))
symlog(gca,'y',0); grid off;
exportgraphics(gca,strcat('Figures/SFS_Alpha4_', num2str(fi), 'xDelta', '.jpeg'),'Resolution',300)
exportgraphics(gca,strcat('Figures/SFS_Alpha4_', num2str(fi), 'xDelta', '.png'),'Resolution',300)

% Figure Pressure Method 2
figure; hold on; box on
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_4_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_4_filt_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_4_closure_FG_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','--','color','k')
plot(y(2,idx:end-idx-1,2)/delta_h,(-1/(v_SFS_FG)*alpha_4_closure_FG_avg(idx:end-idx-1) + alpha_4_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-.','color',[0.4660 0.6740 0.1880])
% plot(y(2,idx:end-idx-1,2)/delta_h,(-1/(v_SFS_FG)*alpha_4_closure_FG_avg(idx:end-idx-1) + alpha_4_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-.','color','red')

xlabel('${y/\delta}$','interpreter','latex')
ylabel(strcat('${\langle (\partial P / \partial t)^\ast \rangle}$', 'Decomposition'),'interpreter','latex')
legend([{'$(\overline{\rho} \breve{c}^2 \nabla \cdot \widetilde{u})^\ast$'},{'$(\overline{{\rho} {c}^2 \nabla \cdot {u}})^\ast$'},...
    {'$({\alpha_4}^{SFS})^\ast$'},{'$(\overline{\rho} \breve{c}^2 \nabla \cdot \widetilde{u} + {\alpha_4}^{SFS})^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.2 1 1])
legend('Location','northeast','box','off','NumColumns', 2)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]); ylim([-1000 1000])
% ylim([10^-2 10^4])
% set(gca,'YTick',10.^(-2:2:4))
symlog(gca,'y',0); grid off;
exportgraphics(gca,strcat('Figures/SFS_Alpha4_', num2str(fi), 'xDelta_FG', '.png'),'Resolution',300)
exportgraphics(gca,strcat('Figures/SFS_Alpha4_', num2str(fi), 'xDelta_FG', '.jpeg'),'Resolution',300)


%% Alpha 5
% alpha_5_unfilt = Alpha_p./(c_v.*Beta_T).*(1./rho.*(Tau_dU - q_div));

alpha_5_SFS_filt   = FilterFields(alpha_5_LES,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);

Alpha_p_filt_filt  = FilterFields(Alpha_p_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
c_v_filt_filt      = FilterFields(c_v_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
Beta_T_filt_filt   = FilterFields(Beta_T_filt,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
Tau_dU_LES_filt    = FilterFields(Tau_dU_LES,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);
q_div_LES_filt     = FilterFields(q_div_LES,delta_filt,delta,fi,x,y,z,dx,dy,dz,Filter_type);

% Strain tensor
divU_filt_favre_favre = du_x_filt_favre_favre./dx + dv_y_filt_favre_favre./dy + dw_z_filt_favre_favre./dz;
alpha_5_closure       = Alpha_p_filt_filt./(c_v_filt_filt.*Beta_T_filt_filt).*(1./rho_filt_filt.*(Tau_dU_LES_filt - q_div_LES_filt)) - alpha_5_SFS_filt;
alpha_5_closure_avg   = Spatial_avg_XZ_var(alpha_5_closure,num_points_x,num_points_y,num_points_z, fi, idx);
alpha_5_filt_avg      = Spatial_avg_XZ_var(alpha_5_filt,num_points_x,num_points_y,num_points_z, fi, idx);

% Constant
v_SFS =  -1/mean((alpha_5_filt_avg(idx:end-idx-1) - alpha_5_LES_avg(idx:end-idx-1))./alpha_5_closure_avg(idx:end-idx-1));

% Figure Pressure
figure; hold on; box on
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_5_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_5_filt_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_5_closure_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','--','color','k')
plot(y(2,idx:end-idx-1,2)/delta_h,(-1/v_SFS*alpha_5_closure_avg(idx:end-idx-1) + alpha_5_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-.','color',[0.4660 0.6740 0.1880])

legend('LES','Closure','Filt','LES + Closure')
xlabel('${y/\delta}$','interpreter','latex')
ylabel(strcat('${\langle (\partial P / \partial t)^\ast \rangle}$', 'Decomposition'),'interpreter','latex')
legend([{'$(\frac{1}{\overline{\rho}} \frac{\breve{\beta_v}}{\breve{c_v} \breve{\beta_T}} ({\breve{\sigma}} : \nabla \otimes \widetilde{\mathbf{u}} - \nabla \cdot \breve{{q}}))^\ast$'},...
    {'$(\overline{\frac{1}{{\rho}} \frac{{\beta_v}}{\breve{c_v} {\beta_T}} ({{\sigma}} : \nabla \otimes {\mathbf{u}} - \nabla \cdot {{q}})})^\ast$'},...
    {'$({\alpha_5}^{SFS})^\ast$'},{'$(\frac{1}{\overline{\rho}} \frac{\breve{\beta_v}}{\breve{c_v} \breve{\beta_T}} ({\breve{\sigma}} : \nabla \otimes \widetilde{\mathbf{u}} - \nabla \cdot \breve{{q}}) + {\alpha_5}^{SFS})^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.2 1 1])
legend('Location','southeast','box','off','NumColumns', 1)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]); ylim([-1000 1000])
% ylim([10^-2 10^4])
% set(gca,'YTick',10.^(-2:2:4))
symlog(gca,'y',0); grid off;
exportgraphics(gca,strcat('Figures/SFS_Alpha5_', num2str(fi), 'xDelta', '.jpeg'),'Resolution',300)
exportgraphics(gca,strcat('Figures/SFS_Alpha5_', num2str(fi), 'xDelta', '.png'),'Resolution',300)



