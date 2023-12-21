%%%%%%%%%%%% Filter TKE assessment plots %%%%%%%%

%% Dataset needs to be loaded for the contourplots!
% Expected to run DNS_filter.m ahead of this script otherwise
% File_name{1} = '3d_high_pressure_turbulent_channel_flow_69100000.h5';
% n_data       = 1;
% File_name_complete = strcat(source_path,File_name{n_data});
% 
% Info      = h5info(File_name_complete);
% Name_var  = {Info.Datasets.Name};
% 
% % Load variables
% for ii = 1:length(Name_var)
%     value = h5read(File_name_complete,strcat('/',Name_var{ii}) );
%     assignin('base',Name_var{ii}, value)
% end

%% Load DataSet
Filter_type = 'CDLF_Box';
name_file_out  = strcat('HP_', 'single_fi_2_10_', Filter_type);
DATA           = readtable(strcat('Data/', name_file_out,'.csv'));

% Load variables
DATA_names = fieldnames(DATA);
for ii = 1:length(DATA_names)
    value = DATA.(DATA_names{ii});
    assignin('base',DATA_names{ii}, value)
end

% Filter width vector
fi           = [1 2 4 6 8 10];

% Dimensions
delta_h   = 100*1E-6; % Channel half height
L_x = 4*pi*delta_h;   % Length in the x direction
L_y = 2*delta_h;      % Length in the y direction
L_z = 4/3*pi*delta_h; % Length in the z direction


%% Plot total TKE
figure; hold on; grid on; box on;
plot(fi,TKE_tot_filt/TKE_tot_filt(1),'o--','linewidth',2,'markersize',5)
xlabel('${{\overline{\Delta}}/\Delta}$','interpreter','latex')
ylabel('${\langle || TKE || \rangle}$','interpreter','latex')

%% Plot filter overshoots
figure; hold on; grid on; box on;

Delta = 4;

FILT = TKE_filt4xDelta(:,64,:);
DNS = TKE(:,64,:);
FILT = FILT(:);
DNS = DNS(:);
plot(FILT(7412:7539),'-','linewidth',2)
plot(DNS(7412:7539),'-','linewidth',2)
xlabel('${N}$','interpreter','latex')
ylabel('${\langle || TKE || \rangle}$','interpreter','latex')
legend([{strcat('${{\overline{\Delta}}/\Delta}$',' = ', num2str(Delta))},{'DNS'}],'interpreter','latex','location','best')
xlim([0 130])
saveas(gca,strcat('Figures/Overshoot_', name_file_out, '_', num2str(Delta), 'xDelta'),'png')
% exportgraphics(gca,strcat('Figures/Overshoot_', name_file_out, '_', num2str(Delta), 'xDelta', '.jpeg'),'Resolution',300)
% exportgraphics(gca,strcat('Figures/Overshoot_', name_file_out, '_', num2str(Delta), 'xDelta', '.png'),'Resolution',300)


%% Plot TKE curves across entire channel y-direction
figure; hold on; grid on; box on;
fi_legend    = cell(1,length(fi));
fi_legend{1} = strcat('${{\overline{\Delta}}/\Delta}$',' =  ', num2str(fi(1)));
plot(y(2:end-1)/delta_h,TKE_y(2:end-1),'-','linewidth',2)
for nn = 2:length(fi)
    QoI = eval(strcat('TKE_y','_filt',num2str(fi(nn)),'xDelta'));
    plot(y(2:end-1)/delta_h,QoI(2:end-1),'-','linewidth',2)
    fi_legend{nn} = strcat('${{\overline{\Delta}}/\Delta}$',' = ', num2str(fi(nn)));
end
ylim([0 max(TKE_y)])
xlabel('${y/\delta}$','interpreter','latex')
ylabel('${\langle TKE \rangle}$','interpreter','latex')
legend(fi_legend,'interpreter','latex','location','best')

%% Plot TKE curves at y/delta position
figure; hold on; grid on; box on
ylim([0 2])

% Set y target
y_delta_target = [0.01 0.25 1 0.75 1.99]; % Define y/delta-position
TKE_y_plot     = zeros(length(y_delta_target),length(fi));

y_legend    = cell(1,length(y_delta_target));
y_legend{1} = strcat('${y/\delta}$',' =  ', num2str(y_delta_target(1)));

for yy = 1:length(y_delta_target)
    [val,idx] = min(abs(y/delta_h-y_delta_target(yy)));
    TKE_y_plot(yy,1) = TKE_y(idx);
    for nn = 2:length(fi)
        TKE_plot          = eval(strcat('TKE_y','_filt',num2str(fi(nn)),'xDelta'));
        TKE_y_plot(yy,nn) = TKE_plot(idx);
    end
    plot(fi,TKE_y_plot(yy,:)./TKE_y_plot(yy,1),'o--','linewidth',2); hold on
    y_legend{yy} = strcat('${y/\delta}$',' = ', num2str(y_delta_target(yy)));
end
xlabel('${{\overline{\Delta}}/\Delta}$','interpreter','latex')
ylabel('${\langle || TKE || \rangle}$','interpreter','latex')
legend(y_legend,'interpreter','latex','location','southwest')
%title(strcat('${\langle || TKE || \rangle}$',' at ', '${y/\delta}$', ' = ', num2str(y_delta_target)),'interpreter','latex')

%% Plot TKE at different y+
figure; hold on; grid on; box on

% Set y target
y_plus_target  = [30 35 40 45 50]; % Define y/delta-position
TKE_y_plot     = zeros(length(y_plus_target),length(fi));

y_legend    = cell(1,length(y_plus_target));
y_legend{1} = strcat('${y/\delta}$',' =  ', num2str(y_plus_target(1)));

% Select bottom (wall = 2:2) or top wall (wall = 1:1)
for wall = 2:2
%     subplot(2,1,wall);
    for yy = 1:length(y_plus_target)
        if wall == 1
            % Top wall
            [val,idx] = min(abs(y_plus_tw - y_plus_target(yy)));
            TKE_y_plot(yy,1) = TKE_tw(idx);
        else
            % Bottom wall
            [val,idx] = min(abs(y_plus_bw - y_plus_target(yy)));
            TKE_y_plot(yy,1) = TKE_bw(idx);
        end
        for nn = 2:length(fi)
            if wall == 1
                % Top wall
                TKE_plot          = eval(strcat('TKE_tw','_filt',num2str(fi(nn)),'xDelta'));
            else
                % Bottom wall
                TKE_plot          = eval(strcat('TKE_bw','_filt',num2str(fi(nn)),'xDelta'));
            end
            TKE_y_plot(yy,nn) = TKE_plot(idx);
        end
        plot(fi,TKE_y_plot(yy,:)./TKE_y_plot(yy,1),'o--','linewidth',2);  hold on; grid on; box on;
        y_legend{yy} = strcat('${y^+}$',' = \thinspace', num2str(y_plus_target(yy)));
    end
    xlabel('${{\overline{\Delta}}/\Delta}$','interpreter','latex')
    ylabel('${\langle || TKE || \rangle}$','interpreter','latex')
    legend(y_legend,'interpreter','latex','location','northeast','box','off')
    ylim([0.3 1]); xlim([1 max(fi)])
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
end

% saveas(gca,strcat('Figures/TKE_vs_filterwidth_y_plus_tw',name_file_out),'png')
exportgraphics(gca,strcat('Figures/TKE_vs_filterwidth_y_plus_bw', name_file_out, '.jpeg'),'Resolution',300)
% exportgraphics(gca,strcat('Figures/TKE_vs_filterwidth_y_plus_tw', name_file_out, '.png'),'Resolution',300)


%% Contourplot streamwise velocity
% DNS
X = x(:,:,64)/delta_h;
Y = y(:,:,64)/delta_h;
Z = u(:,:,64)/u_tau_bw;
f = figure; [c,h]=contourf(X,Y,Z,100);
set(h, 'edgecolor','none'); 
colorbar;
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','u^+$'],'Interpreter','latex','fontsize',10);
xlabel('${{x}/{\delta}}$','interpreter','latex')
ylabel('${{y}/{\delta}}$','interpreter','latex')
xticks(0:2.5:12.5); xticklabels(0:2.5:12.5)
xlim([0 12.5])
yticks(0:0.5:2); yticklabels(0:0.5:2)
ylim([0 2])
set(gca,'linewidth',1)
% set(gca,'fontsize',12)
% ylim([0 1.2])
pbaspect([4 1 1])
% set(gca, 'OuterPosition', [0,0,1,1])
% saveas(gca,strcat('Figures/u_plus_vs_xy',name_file_out,'_',Filter_type,'_DNS'),'png')
exportgraphics(gca,strcat('Figures/u_plus_vs_xy',name_file_out,'_',Filter_type,'_DNS','.jpeg'),'Resolution',300);
% exportgraphics(gca,strcat('Figures/u_plus_vs_xy',name_file_out,'_',Filter_type,'_DNS','.png'),'Resolution',300);

% print(gca, strcat('DNS_Contourplot_u_plus'), 'png', '-r600'); 


% DNS filtered
Z = u_filt2xDelta(:,:,64)/u_tau_bw;
f = figure; [c,h]=contourf(X,Y,Z,100);
% title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
colorbar;
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','u^+$'],'Interpreter','latex','fontsize',10);
xlabel('${{x}/{\delta}}$','interpreter','latex')
ylabel('${{y}/{\delta}}$','interpreter','latex')
xticks(0:2.5:12.5); xticklabels(0:2.5:12.5)
xlim([0 12.5])
yticks(0:0.5:2); yticklabels(0:0.5:2)
ylim([0 2])
set(gca,'linewidth',1)
pbaspect([4 1 1])
set(gca, 'LooseInset', get(gca,'TightInset'))
% saveas(gca,strcat('Figures/u_plus_vs_xy_',name_file_out,'_',Filter_type,'_2xDelta'),'png')
exportgraphics(gca,strcat('Figures/u_plus_vs_xy_',name_file_out,'_',Filter_type,'_2xDelta','.png'),'Resolution',300);
% exportgraphics(gca,strcat('Figures/u_plus_vs_xy_',name_file_out,'_',Filter_type,'_2xDelta','.jpeg'),'Resolution',300);

% DNS filtered
Z = u_filt4xDelta(:,:,64)/u_tau_bw;
f = figure; [c,h]=contourf(X,Y,Z,100);
% title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
caxis([0 15]);colorbar;
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','u^+$'],'Interpreter','latex','fontsize',10);
xlabel('${{x}/{\delta}}$','interpreter','latex')
ylabel('${{y}/{\delta}}$','interpreter','latex')
xticks(0:2.5:12.5); xticklabels(0:2.5:12.5)
xlim([0 12.5])
yticks(0:0.5:2); yticklabels(0:0.5:2)
ylim([0 2])
zlim([0 15]); 
set(gca,'linewidth',1)
pbaspect([4 1 1])
% saveas(gca,strcat('Figures/u_plus_vs_xy_',name_file_out,'_',Filter_type,'_4xDelta'),'png')
exportgraphics(gca,strcat('Figures/u_plus_vs_xy_',name_file_out,'_',Filter_type,'_4xDelta','.jpeg'),'Resolution',300);
% exportgraphics(gca,strcat('Figures/u_plus_vs_xy_',name_file_out,'_',Filter_type,'_4xDelta','.png'),'Resolution',300);


% DNS filtered
Z = u_filt6xDelta(:,:,64)/u_tau_bw;
f = figure; [c,h]=contourf(X,Y,Z,100);
% title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
caxis([0 15]);colorbar;
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','u^+$'],'Interpreter','latex','fontsize',10);
xlabel('${{x}/{\delta}}$','interpreter','latex')
ylabel('${{y}/{\delta}}$','interpreter','latex')
xticks(0:2.5:12.5); xticklabels(0:2.5:12.5)
xlim([0 12.5])
yticks(0:0.5:2); yticklabels(0:0.5:2)
ylim([0 2])
zlim([0 15]); 
set(gca,'linewidth',1)
pbaspect([4 1 1])
% saveas(gca,strcat('Figures/u_plus_vs_xy_',name_file_out,'_',Filter_type,'_6xDelta'),'png')
exportgraphics(gca,strcat('Figures/u_plus_vs_xy_',name_file_out,'_',Filter_type,'_6xDelta','.jpeg'),'Resolution',300);
% exportgraphics(gca,strcat('Figures/u_plus_vs_xy_',name_file_out,'_',Filter_type,'_6xDelta','.png'),'Resolution',300);

%% Contourplot temperature for overleaf
% DNS
Fluid.Substance         = 'N2';
[~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Fluid.Substance);

X = x(:,:,64)/delta_h;
Y = y(:,:,64)/delta_h;
Z = T(:,:,64)/T_c;
f = figure; [c,h]=contourf(X,Y,Z,100);
colormap(redblue);
set(h, 'edgecolor','none'); 
%colorbar;
%cbh = findall(f, 'Type', 'ColorBar');
%cTH = get(cbh,'Title');
%set(cTH,'String',['$','T/T_c$'],'Interpreter','latex','fontsize',10);
%xlabel('${{x}/{\delta}}$','interpreter','latex')
%ylabel('${{y}/{\delta}}$','interpreter','latex')
%xticks(0:2.5:12.5); xticklabels(0:2.5:12.5)
xlim([0 12.5])
%yticks(0.75:0.25:1.5); yticklabels(0.75:0.25:1.5)
ylim([0 2])
% set(gca,'linewidth',1)
set(gca,'xtick',[])
set(gca,'ytick',[])
box off; grid off
% set(gca,'fontsize',12)
pbaspect([4 1 1])
h = gca;
h.YAxis.Visible = 'off';
% set(gca, 'OuterPosition', [0,0,1,1])
% saveas(gca,strcat('Figures/u_plus_vs_xy',name_file_out,'_',Filter_type,'_DNS'),'png')
exportgraphics(gca,strcat('Figures/T_plus_vs_xy',name_file_out,'_',Filter_type,'_DNS','.jpeg'),'Resolution',300);


%% Energy Spectra at different y values

% Spectra 3D
u_inners = u(2:end-1,2:end-1,2:end-1);
v_inners = v(2:end-1,2:end-1,2:end-1);
w_inners = w(2:end-1,2:end-1,2:end-1);

% [k_mag,ke_mag] = Calculate_TKE_spectra_3D(u_inners,v_inners,w_inners,L_x,L_y,L_z,num_points_x-2,num_points_y-2,num_points_z-2);
% figure
% loglog(k_mag,ke_mag,'linewidth',2)
% hold on; grid on; box on;
% loglog(k_mag,ke_mag.^(-5/3),'--')
% legend('DNS','Intertial range')
% xlabel('${{\omega}_{k}}$','interpreter','latex')
% ylabel('${{k_e}({\omega}_{k})}$','interpreter','latex')
% 
% 
% % Bottom wall
% for yy = 1:length(y_vec_norm)
% 
%     % Bottom wall
%     [val,idx] = min(abs(y_plus_bw - y_plus_target(yy)));
% 
% 
%     % For x-direction at first inner Z
%     u_spectra = u(:,idx,2);
%     v_spectra = v(:,idx,2);
%     w_spectra = w(:,idx,2);
%     
% 
% 
% end