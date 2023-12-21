%% FIGURES UNCLOSED TERMS %%

%% Unclosed terms evolution
% figure; hold on; grid on; box on
% subplot(3,2,1)
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_1_avg{1}(idx:end-idx-1)*Norm{1},'linewidth',2)
% hold on
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_1_avg{2}(idx:end-idx-1)*Norm{1},'linewidth',2)
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_1_avg{3}(idx:end-idx-1)*Norm{1},'linewidth',2)
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel('${\langle || \alpha_1 || \rangle}$','interpreter','latex')
% legend([{'x'},{'y'},{'z'}],'interpreter','latex','location','best')
% 
% 
% subplot(3,2,2)
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_2_avg{1}(idx:end-idx-1)*Norm{2},'linewidth',2)
% hold on
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_2_avg{2}(idx:end-idx-1)*Norm{2},'linewidth',2)
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_2_avg{3}(idx:end-idx-1)*Norm{2},'linewidth',2)
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel('${\langle || \alpha_2 || \rangle}$','interpreter','latex')
% legend([{'x'},{'y'},{'z'}],'interpreter','latex','location','best')
% 
% subplot(3,2,3)
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_3_avg{1}(idx:end-idx-1)*Norm{3},'linewidth',2)
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel('${\langle || \alpha_3 || \rangle}$','interpreter','latex')
% 
% subplot(3,2,4)
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_4_avg{1}(idx:end-idx-1)*Norm{4},'linewidth',2)
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel('${\langle || \alpha_4 || \rangle}$','interpreter','latex')
% 
% subplot(3,2,5)
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_5_avg{1}(idx:end-idx-1)*Norm{5},'linewidth',2)
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel('${\langle || \alpha_5 || \rangle}$','interpreter','latex')
% 
% subplot(3,2,6)
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_6_avg{1}(idx:end-idx-1)*Norm{6},'linewidth',2)
% xlabel('${y/\delta}$','interpreter','latex')
% ylabel('${\langle || \alpha_6 || \rangle}$','interpreter','latex')



%% Breakdown
% Momentum X
figure; hold on; box on
% set(gca,'yscale','log')
% plot(y(2,idx:end-idx-1,2)/delta_h,rhou_trans_avg(idx:end-idx-1),'linewidth',2)
plot(y(2,idx:end-idx-1,2)/delta_h,(rhou_conv_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410]); hold on;
plot(y(2,idx:end-idx-1,2)/delta_h,(rhou_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
plot(y(2,idx:end-idx-1,2)/delta_h,(dP_rhou_LES_avg(idx:end-idx-1))*Norm{1},'linewidth',2, 'LineStyle','-','color',[0.8500 0.3250 0.0980])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_1_avg{1}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0 0.4470 0.7410])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_avg{1}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])
xlabel('${y/\delta}$','interpreter','latex')
ylabel(strcat('${\langle (\partial (\overline{\rho} \widetilde{u}) / \partial t)^\ast \rangle}$', ' Decomposition'),'interpreter','latex')
legend([{'$(\nabla \cdot \left( \overline{\rho} \widetilde{\mathbf{u}} \widetilde{\mathbf{u}} \right))^\ast$'},{'$(\nabla \cdot \mathbf{\breve{\sigma}})^\ast$'},{'$(\nabla \overline{P})^\ast$'},{'$\alpha_1^\ast$'},{'$\alpha_2^\ast$'}],'interpreter','latex','location','best')
% legend([{'$(\partial \overline{\rho} \widetilde{u} \widetilde{u} / \partial x)^\ast$'},{'$(\partial \breve{\sigma_x} / \partial x)^\ast$'},{'$(\partial \overline{P} / \partial x)^\ast$'},{'$\alpha_1^\ast$'},{'$\alpha_2^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.2 1 1])
legend('Location','northwest','box','off','NumColumns', 2)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]); ylim([-0.1 0.1])
symlog(gca,'y',-4); grid off;
% ylim([-10^0 10^0])
% set(gca,'YTick',10.^(-5:0))
% set(gca,'YMinorTick','on')
% saveas(gca,strcat('Figures/Unclosed_Momentum_X_', num2str(fi), 'xDelta'),'png')
% exportgraphics(gca,strcat('Figures/Unclosed_Momentum_X_', num2str(fi), 'xDelta', '.png'),'Resolution',300)
exportgraphics(gca,strcat('Figures/Unclosed_Momentum_X_', num2str(fi), 'xDelta', '.jpeg'),'Resolution',300)

% Momentum Y
figure; hold on; box on
% set(gca,'yscale','log')
% plot(y(2,idx:end-idx-1,2)/delta_h,rhou_trans_avg(idx:end-idx-1),'linewidth',2)
plot(y(2,idx:end-idx-1,2)/delta_h,(rhov_conv_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410]); hold on;
plot(y(2,idx:end-idx-1,2)/delta_h,(rhov_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
plot(y(2,idx:end-idx-1,2)/delta_h,(dP_rhov_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.8500 0.3250 0.0980])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_1_avg{2}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0 0.4470 0.7410])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_avg{2}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])
xlabel('${y/\delta}$','interpreter','latex')
ylabel(strcat('${\langle (\partial (\overline{\rho} \widetilde{v}) / \partial t)^\ast \rangle}$', ' Decomposition'),'interpreter','latex')
legend([{'$(\nabla \cdot \left( \overline{\rho} \widetilde{\mathbf{u}} \widetilde{\mathbf{u}} \right))^\ast$'},{'$(\nabla \cdot \mathbf{\breve{\sigma}})^\ast$'},{'$(\nabla \overline{P})^\ast$'},{'$\alpha_1^\ast$'},{'$\alpha_2^\ast$'}],'interpreter','latex','location','best')
% legend([{'$(\partial \overline{\rho} \widetilde{u} \widetilde{v} / \partial y)^\ast$'},{'$(\partial \breve{\sigma_y} / \partial y)^\ast$'},{'$(\partial \overline{P} / \partial y)^\ast$'},{'$\alpha_1^\ast$'},{'$\alpha_2^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.2 1 1])
legend('Location','northwest','box','off','NumColumns', 2)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]);  ylim([-0.1 0.1])
% ylim([10^-2 10^0])
symlog(gca,'y',-4); grid off;

% set(gca,'YTick',10.^(-8:2:0))
% saveas(gca,strcat('Figures/Unclosed_Momentum_Y_', num2str(fi), 'xDelta'),'png')
% exportgraphics(gca,strcat('Figures/Unclosed_Momentum_Y_', num2str(fi), 'xDelta', '.png'),'Resolution',300)
exportgraphics(gca,strcat('Figures/Unclosed_Momentum_Y_', num2str(fi), 'xDelta', '.jpeg'),'Resolution',300)


% Momentum Z
figure; hold on; box on
% set(gca,'yscale','log')
% plot(y(2,idx:end-idx-1,2)/delta_h,rhou_trans_avg(idx:end-idx-1),'linewidth',2)
plot(y(2,idx:end-idx-1,2)/delta_h,(rhow_conv_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410]); hold on;
plot(y(2,idx:end-idx-1,2)/delta_h,(rhow_visc_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
plot(y(2,idx:end-idx-1,2)/delta_h,(dP_rhow_LES_avg(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','-','color',[0.8500 0.3250 0.0980])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_1_avg{3}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0 0.4470 0.7410])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_2_avg{3}(idx:end-idx-1))*Norm{1},'LineWidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])
xlabel('${y/\delta}$','interpreter','latex')
ylabel(strcat('${\langle (\partial (\overline{\rho} \widetilde{w}) / \partial t)^\ast \rangle}$', ' Decomposition'),'interpreter','latex')
legend([{'$(\nabla \cdot \left( \overline{\rho} \widetilde{\mathbf{u}} \widetilde{\mathbf{u}} \right))^\ast$'},{'$(\nabla \cdot \mathbf{\breve{\sigma}})^\ast$'},{'$(\nabla \overline{P})^\ast$'},{'$\alpha_1^\ast$'},{'$\alpha_2^\ast$'}],'interpreter','latex','location','best')
% legend([{'$(\partial \overline{\rho} \widetilde{u} \widetilde{w} / \partial z)^\ast$'},{'$(\partial \breve{\sigma_z} / \partial z)^\ast$'},{'$(\partial \overline{P} / \partial z)^\ast$'},{'$\alpha_1^\ast$'},{'$\alpha_2^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.2 1 1])
legend('Location','north','box','off','NumColumns', 2)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]); ylim([-0.01 0.01]) %ylim([10^-8 10^0])
% set(gca,'YTick',10.^(-8:2:0))
symlog(gca,'y',-4); grid off;
% saveas(gca,strcat('Figures/Unclosed_Momentum_Z_', num2str(fi), 'xDelta'),'png')
% exportgraphics(gca,strcat('Figures/Unclosed_Momentum_Z_', num2str(fi), 'xDelta', '.png'),'Resolution',300)
exportgraphics(gca,strcat('Figures/Unclosed_Momentum_Z_', num2str(fi), 'xDelta', '.jpeg'),'Resolution',300)


% Pressure
figure; hold on; box on
% set(gca,'yscale','log')
% plot(y(2,idx:end-idx-1,2)/delta_h,P_trans_avg(idx:end-idx-1),'linewidth',2)
plot(y(2,idx:end-idx-1,2)/delta_h,(u_grad_P_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-','color',[0 0.4470 0.7410]); hold on;
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_4_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_5_LES_avg(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','-','color',[0.8500 0.3250 0.0980])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_3_avg{1}(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','--','color',[0 0.4470 0.7410])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_4_avg{1}(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','--','color',[0.4660 0.6740 0.1880])
plot(y(2,idx:end-idx-1,2)/delta_h,(alpha_5_avg{1}(idx:end-idx-1))*Norm{3},'linewidth',2, 'LineStyle','--','color',[0.8500 0.3250 0.0980])
xlabel('${y/\delta}$','interpreter','latex')
ylabel(strcat('${\langle (\partial P / \partial t)^\ast \rangle}$', 'Decomposition'),'interpreter','latex')
legend([{'$(\widetilde{u} \nabla \overline{P})^\ast$'},{'$(\overline{\rho} \breve{c}^2 \nabla \cdot \widetilde{u})^\ast$'},{'$(\frac{1}{\overline{\rho}} \frac{\breve{\beta_v}}{\breve{c_v} \breve{\beta_T}} ({\breve{\sigma}} : \nabla \otimes \widetilde{\mathbf{u}} - \nabla \cdot \breve{{q}}))^\ast$'},{'$\alpha_3^\ast$'},{'$\alpha_4^\ast$'},{'$\alpha_5^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.2 1 1])
legend('Location','northeast','box','off','NumColumns', 2)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]); ylim([-1000 1000])
% ylim([10^-2 10^4])
% set(gca,'YTick',10.^(-2:2:4))
symlog(gca,'y',0); grid off;

% saveas(gca,strcat('Figures/Unclosed_Pressure_', num2str(fi), 'xDelta'),'png')
% exportgraphics(gca,strcat('Figures/Unclosed_Pressure_', num2str(fi), 'xDelta', '.png'),'Resolution',300)
exportgraphics(gca,strcat('Figures/Unclosed_Pressure_', num2str(fi), 'xDelta', '.jpeg'),'Resolution',300)


% EOS
figure; hold on; box on
% set(gca,'yscale','log')
semilogy(y(2,idx:end-idx-1,2)/delta_h,(P_LES_avg(idx:end-idx-1))*Norm{6},'linewidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
hold on
semilogy(y(2,idx:end-idx-1,2)/delta_h,(P_filt_avg(idx:end-idx-1))*Norm{6},'linewidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_6_avg{1}(idx:end-idx-1),'linewidth',2)
xlabel('${y/\delta}$','interpreter','latex')
ylabel('${\langle P^\ast \rangle}$','interpreter','latex')
legend([{'$P(\overline{\rho},\breve{T})^\ast$'},{'$\overline{P(\rho,T)}^\ast$'}],'interpreter','latex','location','best')
pbaspect([1.8 1 1])
legend('Location','northwest','box','off','NumColumns', 3)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 2]); ylim([6200 6500])
% axis tight
% saveas(gca,strcat('Figures/Unclosed_EOS_', num2str(fi), 'xDelta'),'png')
% exportgraphics(gca,strcat('Figures/Unclosed_EOS_', num2str(fi), 'xDelta', '.png'),'Resolution',300)
exportgraphics(gca,strcat('Figures/Unclosed_EOS_', num2str(fi), 'xDelta', '.jpeg'),'Resolution',300)



%% Obtain data for y-positions
% Max stremwise velocity position
[avg_u_xz]     = Spatial_avg_XZ_var(avg_u,num_points_x,num_points_y,num_points_z, 1);
u_max          = max(avg_u_xz);
[val,idx_umax] = min(abs(avg_u_xz-u_max));
y_umax = y(1,idx_umax,1)/delta_h;
y_umax = 1;

y_delta_target = [0.2 y_umax 1.8]; % Define y/delta-position

for yy = 1:length(y_delta_target)
    [val,idx_y(yy)] = min(abs(y(1,:,1)/delta_h-y_delta_target(yy)));
end


disp("Momentum X unclosed values ...")
disp("rhou_conv_LES_avg = " + num2str(rhou_conv_LES_avg(idx_y)*Norm{1},'%e'))
disp("rhou_visc_LES_avg = " + num2str(rhou_visc_LES_avg(idx_y)*Norm{1},'%e'))
disp("dP_rhou_LES_avg = " + num2str(dP_rhou_LES_avg(idx_y)*Norm{1},'%e'))
disp("alpha_1_avg{1} = " + num2str(alpha_1_avg{1}(idx_y)*Norm{1},'%e'))
disp("alpha_2_avg{1} = " + num2str(alpha_2_avg{1}(idx_y)*Norm{1},'%e'))
disp("ratio alpha_1 = " + num2str(alpha_1_avg{1}(idx_y)./(-rhou_conv_LES_avg(idx_y) + rhou_visc_LES_avg(idx_y) - dP_rhou_LES_avg(idx_y)),'%e'))
disp("ratio alpha_2 = " + num2str(alpha_2_avg{1}(idx_y)./(-rhou_conv_LES_avg(idx_y) + rhou_visc_LES_avg(idx_y) - dP_rhou_LES_avg(idx_y)),'%e'))

disp(' ')

disp("Momentum Y unclosed values ...")
disp("rhou_conv_LES_avg = " + num2str(rhov_conv_LES_avg(idx_y)*Norm{1},'%e'))
disp("rhou_visc_LES_avg = " + num2str(rhov_visc_LES_avg(idx_y)*Norm{1},'%e'))
disp("dP_rhou_LES_avg = " + num2str(dP_rhov_LES_avg(idx_y)*Norm{1},'%e'))
disp("alpha_1_avg{2} = " + num2str(alpha_1_avg{2}(idx_y)*Norm{1},'%e'))
disp("alpha_2_avg{2} = " + num2str(alpha_2_avg{2}(idx_y)*Norm{1},'%e'))
disp("ratio alpha_1 = " + num2str(alpha_1_avg{2}(idx_y)./(-rhov_conv_LES_avg(idx_y) + rhov_visc_LES_avg(idx_y) - dP_rhov_LES_avg(idx_y)),'%e'))
disp("ratio alpha_2 = " + num2str(alpha_2_avg{2}(idx_y)./(-rhov_conv_LES_avg(idx_y) + rhov_visc_LES_avg(idx_y) - dP_rhov_LES_avg(idx_y)),'%e'))

disp(' ')

disp("Momentum Z unclosed values ...")
disp("rhou_conv_LES_avg = " + num2str(rhow_conv_LES_avg(idx_y)*Norm{1},'%e'))
disp("rhou_visc_LES_avg = " + num2str(rhow_visc_LES_avg(idx_y)*Norm{1},'%e'))
disp("dP_rhou_LES_avg = " + num2str(dP_rhow_LES_avg(idx_y)*Norm{1},'%e'))
disp("alpha_1_avg{3} = " + num2str(alpha_1_avg{3}(idx_y)*Norm{1},'%e'))
disp("alpha_2_avg{3} = " + num2str(alpha_2_avg{3}(idx_y)*Norm{1},'%e'))
disp("ratio alpha_1 = " + num2str(alpha_1_avg{3}(idx_y)./(-rhow_conv_LES_avg(idx_y) + rhow_visc_LES_avg(idx_y) - dP_rhow_LES_avg(idx_y)),'%e'))
disp("ratio alpha_2 = " + num2str(alpha_2_avg{3}(idx_y)./(-rhow_conv_LES_avg(idx_y) + rhow_visc_LES_avg(idx_y) - dP_rhow_LES_avg(idx_y)),'%e'))

disp(' ')

disp("Pressure unclosed values ...")
F_tot = (-u_grad_P_LES_avg(idx) - alpha_4_LES_avg(idx_y) + alpha_5_LES_avg(idx_y));
disp("F_tot = " + num2str(F_tot*Norm{3}),'%e')
disp("alpha_3_avg{1} = " + num2str(alpha_3_avg{1}(idx_y)*Norm{3},'%e'))
disp("alpha_4_avg{1} = " + num2str(alpha_4_avg{1}(idx_y)*Norm{3},'%e'))
disp("alpha_5_avg{1} = " + num2str(alpha_5_avg{1}(idx_y)*Norm{3},'%e'))
disp("ratio alpha_3 = " + num2str(alpha_3_avg{1}(idx_y)./(F_tot),'%e'))
disp("ratio alpha_4 = " + num2str(alpha_4_avg{1}(idx_y)./(F_tot),'%e'))
disp("ratio alpha_5 = " + num2str(alpha_5_avg{1}(idx_y)./(F_tot),'%e'))

disp(' ')

disp("EOS unclosed values ...")
disp("P_LES_avg = " + num2str(P_LES_avg(idx_y)*Norm{6},'%e'))
disp("P_filt_avg = " + num2str(P_filt_avg(idx_y)*Norm{6},'%e'))
disp("alpha_6_avg{1} = " + num2str(alpha_6_avg{1}(idx_y)*Norm{6},'%e'))
disp("ratio alpha_6 = " + num2str(alpha_6_avg{1}(idx_y)./(P_LES_avg(idx_y)),'%e'))

disp(' ')
