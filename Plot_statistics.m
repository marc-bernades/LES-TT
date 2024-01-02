%% Plot statistics

% Load results from data set and average

%% Load DataSet
Filter_type = 'CDLF_Box';
name_file_out  = strcat('Data/', 'HP_fi_2_10_single', Filter_type);
% name_file_out  = strcat('Data/', 'LP_', Filter_type);

DATA           = readtable(strcat(name_file_out,'.csv'));

% Dimensions
delta_h   = 100*1E-6; % Channel half height
L_x = 4*pi*delta_h;   % Length in the x direction
L_y = 2*delta_h;      % Length in the y direction
L_z = 4/3*pi*delta_h; % Length in the z direction

% Load HP
% num_points_x = unique(DATA.X);
% num_points_y = unique(DATA.Y);
% num_points_z = unique(DATA.Z);
% u         = reshape(DATA.u,[num_points_x,num_points_y,num_points_z]);


%% 1st order statistics
% Maximum velocity Bottom wall
max_y_index_bw = -1;
max_value_bw      = -1.0;
for p = 1:length( DATA.y_plus_bw )
    value = DATA.u_plus_bw(p);
    if( value > max_value_bw )
        max_y_index_bw = p;
        max_value_bw   = value;
    else
        disp("y+_bw position at max velocity = " + num2str(DATA.y_plus_bw(max_y_index_bw)) )
        break
    end
end

% Maximum velocity Top wall
max_y_index_tw = -1;
max_value_tw      = -1.0;
for p = 1:length( DATA.y_plus_tw )
    value = DATA.u_plus_tw(p);
    if( value > max_value_tw )
        max_y_index_tw = p;
        max_value_tw   = value;
    else
        disp("y+_tw position at max velocity = " + num2str(DATA.y_plus_tw(max_y_index_tw)) )
        break
    end
end


% 1st order statistics (velocity) high-pressure (DNS, filter 2x and 4x)
figure; hold on; grid on; box on;
xticklabels([]); yticklabels([])
h = axes;
set(h,'xscale','log')
semilogx(DATA.y_plus_bw(2:max_y_index_bw),DATA.u_plus_bw(2:max_y_index_bw),'o','LineWidth',1.5,'MarkerSize',8,'color',[0 0.4470 0.7410]); hold on;grid on; box on;
semilogx(DATA.y_plus_bw(2:max_y_index_bw),DATA.u_plus_bw_filt2xDelta(2:max_y_index_bw),'s','LineWidth',1.5,'MarkerSize',8,'color',[0.6350 0.0780 0.1840])
semilogx(DATA.y_plus_bw(2:max_y_index_bw),DATA.u_plus_bw_filt4xDelta(2:max_y_index_bw),'x','LineWidth',1.5,'MarkerSize',8,'color',[0.4660 0.6740 0.1880])
% ylim([0 15])
xticks([10^-1 10^0 10^1 10^2])
xlabel('${y^+}$','interpreter','latex')
ylabel('${u^+}$','interpreter','latex')
legend([{'DNS'},{strcat('${{\overline{\Delta}}/\Delta}$',' = ','$\thinspace$', num2str(2))},{strcat('${{\overline{\Delta}}/\Delta}$',' = ', '$\thinspace$',num2str(4))}],'interpreter','latex','location','best')
legend('Location','best','box','off')
set(gca,'linewidth',2)
set(gca,'fontsize',16)
exportgraphics(gca,strcat('Figures/uplus_vs_yplus_bw_', name_file_out(6:end), '.jpeg'),'Resolution',300)

% 1st order statistics (velocity) high-pressure (DNS, filter 2x and 4x)
figure; hold on; grid on; box on;
xticklabels([]); yticklabels([])
h = axes;
set(h,'xscale','log')
semilogx(DATA.y_plus_tw(2:max_y_index_tw),DATA.u_plus_tw(2:max_y_index_tw),'o','LineWidth',1.5,'MarkerSize',8,'color',[0 0.4470 0.7410]); hold on;grid on; box on;
semilogx(DATA.y_plus_tw(2:max_y_index_tw),DATA.u_plus_tw_filt2xDelta(2:max_y_index_tw),'s','LineWidth',1.5,'MarkerSize',8,'color',[0.6350 0.0780 0.1840])
semilogx(DATA.y_plus_tw(2:max_y_index_tw),DATA.u_plus_tw_filt4xDelta(2:max_y_index_tw),'x','LineWidth',1.5,'MarkerSize',8,'color',[0.4660 0.6740 0.1880])
% ylim([0 20])
xticks([10^-1 10^0 10^1 10^2])
xlabel('${y^+}$','interpreter','latex','fontsize',16)
ylabel('${u^+}$','interpreter','latex','fontsize',16)
legend([{'DNS'},{strcat('${{\overline{\Delta}}/\Delta}$',' = ', '$\thinspace$', num2str(2))},{strcat('${{\overline{\Delta}}/\Delta}$',' = ','$\thinspace$', num2str(4))}],'interpreter','latex','location','best')
legend('Location','best','box','off','fontsize',16)
set(gca,'linewidth',2)
set(gca,'fontsize',16)
exportgraphics(gca,strcat('Figures/uplus_vs_yplus_tw_', name_file_out(6:end), '.jpeg'),'Resolution',300)


%% u/u_b vs. y/delta
figure; hold on; grid on; box on;
plot(DATA.y/delta_h,DATA.avg_u_xz./DATA.u_b,'o','LineWidth',1.5,'MarkerSize',6,'color',[0 0.4470 0.7410]);
plot(DATA.y/delta_h,DATA.avg_u_xz_filt2xDelta./DATA.u_b_filt2xDelta,'x','LineWidth',1.5,'MarkerSize',6,'color',[0.6350 0.0780 0.1840]);
plot(DATA.y/delta_h,DATA.avg_u_xz_filt4xDelta./DATA.u_b_filt4xDelta,'+','LineWidth',1.5,'MarkerSize',6,'color',[0.4660 0.6740 0.1880]);
xlim([0 2])
xlabel('${y/\delta}$','interpreter','latex')
ylabel('${\langle u/{u_b} \rangle }$','interpreter','latex')
legend([{'DNS'},{strcat('${{\overline{\Delta}}/\Delta}$',' = ', num2str(2))},{strcat('${{\overline{\Delta}}/\Delta}$',' = ', num2str(4))}],'interpreter','latex','location','best')
legend('Location','best','box','off')
set(gca,'linewidth',2)
set(gca,'fontsize',16)
% saveas(gca,strcat('Figures/U_vs_y_',name_file_out),'png')
exportgraphics(gca,strcat('Figures/U_vs_y_',name_file_out(6:end), '.jpeg'),'Resolution',300)
% exportgraphics(gca,strcat('Figures/U_vs_y_',name_file_out, '.png'),'Resolution',300)


%% u_rms/u_b vs y/delta
figure; hold on; grid on; box on;
plot(DATA.y/delta_h,sqrt(DATA.avg_R_favre_uu_bw)./DATA.u_b,'o','LineWidth',1.5,'MarkerSize',5,'color',[0 0.4470 0.7410]);
plot(DATA.y/delta_h,sqrt(DATA.avg_R_favre_uu_bw_filt2xDelta)./DATA.u_b_filt2xDelta,'o','LineWidth',1.5,'MarkerSize',5,'color',[0.6350 0.0780 0.1840]);
plot(DATA.y/delta_h,sqrt(DATA.avg_R_favre_uu_bw_filt4xDelta)./DATA.u_b_filt4xDelta,'o','LineWidth',1.5,'MarkerSize',5,'color',[0.4660 0.6740 0.1880]);
plot(DATA.y/delta_h,sqrt(DATA.avg_R_favre_vv_bw)./DATA.u_b,'+','LineWidth',1.5,'MarkerSize',5,'color',[0 0.4470 0.7410]);
plot(DATA.y/delta_h,sqrt(DATA.avg_R_favre_vv_bw_filt2xDelta)./DATA.u_b_filt2xDelta,'+','LineWidth',1.5,'MarkerSize',5,'color',[0.6350 0.0780 0.1840]);
plot(DATA.y/delta_h,sqrt(DATA.avg_R_favre_vv_bw_filt4xDelta)./DATA.u_b_filt4xDelta,'+','LineWidth',1.5,'MarkerSize',5,'color',[0.4660 0.6740 0.1880]);
plot(DATA.y/delta_h,sqrt(DATA.avg_R_favre_ww_bw)./DATA.u_b,'s','LineWidth',1.5,'MarkerSize',5,'color',[0 0.4470 0.7410]);
plot(DATA.y/delta_h,sqrt(DATA.avg_R_favre_ww_bw_filt2xDelta)./DATA.u_b_filt2xDelta,'s','LineWidth',1.5,'MarkerSize',5,'color',[0.6350 0.0780 0.1840]);
plot(DATA.y/delta_h,sqrt(DATA.avg_R_favre_ww_bw_filt4xDelta)./DATA.u_b_filt4xDelta,'s','LineWidth',1.5,'MarkerSize',5,'color',[0.4660 0.6740 0.1880]);

xlim([0 2]); ylim([0 0.25])
xlabel('${y/\delta}$','interpreter','latex')
ylabel('${\langle u_{rms}/{u_b} \rangle, \langle v_{rms}/{u_b} \rangle, \langle w_{rms}/{u_b} \rangle}$','interpreter','latex')
legend([{'${u_{rms}}  \thinspace {DNS}$'},{strcat('${u_{rms}}', ' \thinspace {{\overline{\Delta}}/\Delta}$',' = ', num2str(2))},{strcat('${u_{rms}}', ' \thinspace {{\overline{\Delta}}/\Delta}$',' = ', num2str(4))}, ...
    {'${v_{rms}} \thinspace {DNS}$'},{strcat('${v_{rms}}', ' \thinspace {{\overline{\Delta}}/\Delta}$',' = ', num2str(2))},{strcat('${v_{rms}}', ' \thinspace {{\overline{\Delta}}/\Delta}$',' = ', num2str(4))}, ...
    {'${w_{rms}}  \thinspace {DNS}$'},{strcat('${w_{rms}}', ' \thinspace {{\overline{\Delta}}/\Delta}$',' = ', num2str(2))},{strcat('${w_{rms}}', ' \thinspace {{\overline{\Delta}}/\Delta}$',' = ', num2str(4))}],'interpreter','latex','location','best')
legend('Location','best','box','off', 'NumColumns', 3)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
pbaspect([2.0 1.5 1])
% set(gca, 'OuterPosition', [0,0,1,1])
% saveas(gca,strcat('Figures/RMS_',name_file_out),'png')
exportgraphics(gca,strcat('Figures/RMS_',name_file_out(6:end), '.jpeg'),'Resolution',300)
% exportgraphics(gca,strcat('Figures/RMS_',name_file_out, '.png'),'Resolution',300)


