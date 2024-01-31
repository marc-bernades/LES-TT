function EigenDecomposition_SFS_Model(x,y,z,idx_filt, Tau_xx,Tau_xy,Tau_xz,Tau_yy,Tau_yz,Tau_zz,Tau_kk_u_SFS, Tau_kk_v_SFS, Tau_kk_w_SFS,Tau_kk_SFS,u_b,delta_h, fi, SFS_model)

%% Eigen Decomposition
% close all; clc

% Initialization
num_points_x = length(x(:,1,1));
num_points_y = length(x(1,:,1));
num_points_z = length(x(1,1,:));

% Barycentric map corners
x1c = [0, 0];
x2c = [1, 0];
x3c = [1/2 sqrt(3)/2];

% Y - positions for barycentric map
y_outer_vec = [0.2 1 1.8]*delta_h;

for yy = 1:length(y_outer_vec)

    [~,idx_y]    = min(abs(y(1,:,1)-y_outer_vec(yy)));

    bar_x_map  = [];
    bar_y_map  = [];
    bar_kk_map = [];

    radius_map_1 = [];  radius_map_2 = [];  radius_map_3 = [];
    phi_map_1    = [];  phi_map_2    = [];  phi_map_3   = [];
    theta_map_1  = [];  theta_map_2  = [];  theta_map_3 = [];


    % Decomposition point by point
    for jj = idx_y:idx_y
        for ii = idx_filt:num_points_x-(idx_filt-1) % Only filtered points from y_outer
            for kk = idx_filt:num_points_z-(idx_filt-1)

                tau_filt_ij = [Tau_xx(ii,jj,kk), Tau_xy(ii,jj,kk), Tau_xz(ii,jj,kk);
                    Tau_xy(ii,jj,kk), Tau_yy(ii,jj,kk), Tau_yz(ii,jj,kk);
                    Tau_xz(ii,jj,kk), Tau_yz(ii,jj,kk), Tau_zz(ii,jj,kk)];


                % From Trace Model
                tau_filt_ij(1,1) = tau_filt_ij(1,1) + Tau_kk_u_SFS(ii,jj,kk);
                tau_filt_ij(2,2) = tau_filt_ij(2,2) + Tau_kk_v_SFS(ii,jj,kk);
                tau_filt_ij(3,3) = tau_filt_ij(3,3) + Tau_kk_w_SFS(ii,jj,kk);

                tau_filt_kk = tau_filt_ij(1,1) + tau_filt_ij(2,2) + tau_filt_ij(3,3);

%                 tau_filt_kk = Tau_kk_SFS(ii,jj,kk);

                % Remove back trace
                tau_filt_ij(1,1) = tau_filt_ij(1,1) - 1/3*tau_filt_kk;
                tau_filt_ij(2,2) = tau_filt_ij(2,2) - 1/3*tau_filt_kk;
                tau_filt_ij(3,3) = tau_filt_ij(3,3) - 1/3*tau_filt_kk;

                % Normalize
                tau_filt_ij_norm = tau_filt_ij./tau_filt_kk;

                [tau_eigenvec_raw, tau_D_eigenval] = eig(tau_filt_ij_norm);
                tau_eigenval_raw                   = eig(tau_filt_ij_norm);
                [tau_eigenval, idx_sort]           = sort(tau_eigenval_raw,'descend'); % First greatest eigenvalue
                tau_eigenvec                       = tau_eigenvec_raw(:, idx_sort);    % Sort based on eigenvalue

                % Prepare barycentric map coordinates
                bar_xy_map = x1c.*(tau_eigenval(1) - tau_eigenval(2)) + 2*x2c.*(tau_eigenval(2) - tau_eigenval(3)) + x3c.*(3*tau_eigenval(3) + 1);

                bar_x_map  = [bar_x_map, bar_xy_map(1)];
                bar_y_map  = [bar_y_map, bar_xy_map(2)];
                bar_kk_map = [bar_kk_map, tau_filt_kk];

                % Prepare eigenvector radius and angles

                % Radius
                radius{1} = sqrt(tau_eigenvec(1,1).^2 + tau_eigenvec(2,1).^2 + tau_eigenvec(3,1).^2);
                radius{2} = sqrt(tau_eigenvec(1,2).^2 + tau_eigenvec(2,2).^2 + tau_eigenvec(3,2).^2);
                radius{3} = sqrt(tau_eigenvec(1,3).^2 + tau_eigenvec(2,3).^2 + tau_eigenvec(3,3).^2);

                radius_map_1 = [radius_map_1, radius{1}];
                radius_map_2 = [radius_map_2, radius{2}];
                radius_map_3 = [radius_map_3, radius{3}];

                % Polar angle (0 - pi)
                phi{1} = acos(tau_eigenvec(3,1)/radius{1});
                phi{2} = acos(tau_eigenvec(3,2)/radius{2});
                phi{3} = acos(tau_eigenvec(3,3)/radius{3});

                phi_map_1 = [phi_map_1, phi{1}];
                phi_map_2 = [phi_map_2, phi{2}];
                phi_map_3 = [phi_map_3, phi{3}];

                % Theta angle (azimuth: [-pi,pi])
                theta{1} = atan(tau_eigenvec(2,1)/tau_eigenvec(1,1));
                theta{2} = atan(tau_eigenvec(2,2)/tau_eigenvec(1,2));
                theta{3} = atan(tau_eigenvec(2,3)/tau_eigenvec(1,3));

                theta_map_1 = [theta_map_1, theta{1}];
                theta_map_2 = [theta_map_2, theta{2}];
                theta_map_3 = [theta_map_3, theta{3}];


            end

        end
    end



    n_bins_pdf          = 2*ceil(1 + log2(length(bar_kk_map)));


    %% PDF radius (eigenvector distance)
    radius_map_1_norm{yy} = radius_map_1./delta_h;
    radius_map_2_norm{yy} = radius_map_2./delta_h;
    radius_map_3_norm{yy} = radius_map_3./delta_h;

    n_bins_pdf_radius_1   = 2*ceil(1 + log2(length(radius_map_1)));
    n_bins_pdf_radius_2   = 2*ceil(1 + log2(length(radius_map_2)));
    n_bins_pdf_radius_3   = 2*ceil(1 + log2(length(radius_map_3)));

    %% PDF on barycentric map (eigen values)
    % Prepare histogram
    n_bins = [50 50];

    [N,Xedges,Yedges] = histcounts2(bar_x_map,bar_y_map,n_bins,'Xbinlimits',[0 1],'Ybinlimits',[0 1],'Normalization','pdf');

    N = N + 1E-12;

    % Remove outlier N
%     if max(max(N)) > 1E3
%         [~, idx_max_j] = max(max(N));
%         [~,idx_max_i] = max(N(:,idx_max_j));
%         N(idx_max_i, idx_max_j) = NaN;
%     end

    x_center = 0.5*(Xedges(2:end) + Xedges(1:end-1));
    y_center = 0.5*(Yedges(2:end) + Yedges(1:end-1));

    % Extend center point to exact to barycentric limit
    x_center = [0 x_center 1];
    y_center = [0 y_center x3c(2)];
    N = [zeros(1,length(N(1,:))); N; zeros(1,length(N(1,:)))];
    N = [zeros(length(N(:,1)),1) N zeros(length(N(:,1)),1)];

    f = figure;
    h = pcolor(x_center,y_center,N');
    set(h, 'EdgeColor', 'none');
    colorbar; colormap((jet)) % flipud
    cbh = findall(f, 'Type', 'ColorBar');
    cTH = get(cbh,'Title');
    set(cTH,'String',['$','PDF','$'],'Interpreter','latex','fontsize',16);
    set(gca,'linewidth',1)
    hold on
    %     x_lim_white = 1/tan(60*pi/180); % y = 1
    %     fill([0 -0.1 -0.1 x_lim_white], [0 0 1 1], 'white', 'LineStyle','none')
    %     fill([1-x_lim_white 1.1 1.1 1.0], [1 1 0 0], 'white', 'LineStyle','none')

    % Remove colors outside barycentric region
    fill([0 -0.1 -0.1 0.5], [0 0 x3c(2) x3c(2)], 'white', 'LineStyle','none')
    fill([0.5 1.1 1.1 1.0], [x3c(2) x3c(2) 0 0], 'white', 'LineStyle','none')
    fill([0 0 1 1], [x3c(2) x3c(2)+0.2 x3c(2)+0.2 x3c(2)], 'white', 'LineStyle','none')
    xlim([0 1])
    ylim([0 x3c(2)])
    axis off

    % Add barycentric vertex junctions
    hold on
    plot([x2c(1) x1c(1)],[x2c(2) x1c(2)],'k','LineWidth',1)
    plot([x3c(1) x1c(1)],[x3c(2) x1c(2)],'k','LineWidth',1)
    plot([x3c(1) x2c(1)],[x3c(2) x2c(2)],'k','LineWidth',1)

    box off; axis off;
    text(x1c(1)-0.05, x1c(2)-0.02, '${x_1}_c$', 'interpreter', 'latex','fontsize',16)
    text(x2c(1), x2c(2)-0.02, '${x_2}_c$', 'interpreter', 'latex','fontsize',16)
    text(x3c(1)-0.02, x3c(2)+0.02, '${x_3}_c$', 'interpreter', 'latex','fontsize',16)

%     saveas(gca,strcat('Figures/PDF_Barycentric_map_DNS_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model),'epsc')
%     saveas(gca,strcat('Figures/PDF_Barycentric_map_DNS_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model),'png')
%     saveas(gca,strcat('Figures/PDF_Barycentric_map_DNS_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model),'eps')
%     exportgraphics(gca,strcat('Figures/PDF_Barycentric_map_DNS_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model, '.png'),'Resolution',300)
    exportgraphics(gca,strcat('Figures/PDF_Barycentric_map_DNS_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model, '.jpeg'),'Resolution',300)

    %% PDF on Polar map (eigen vectors)
    % Prepare histogram
    figure;
    polarhistogram(phi_map_1,n_bins_pdf, 'DisplayStyle','stairs','LineWidth',2, 'LineStyle','-','EdgeColor',[0 0.4470 0.7410]); %'Normalization', 'pdf'
    hold on
    polarhistogram(phi_map_2,n_bins_pdf,'DisplayStyle','stairs', 'LineWidth',2, 'LineStyle','-','EdgeColor',[0.4660 0.6740 0.1880]); %'Normalization', 'pdf'
    polarhistogram(phi_map_3,n_bins_pdf,'DisplayStyle','stairs', 'LineWidth',2, 'LineStyle','-','EdgeColor',[0.8500 0.3250 0.0980]); %'Normalization', 'pdf'
    thetalim([0 180])
    legend([{'$\phi_1$'},{'$\phi_2$'},{'$\phi_3$'}],'interpreter','latex','location','southwest','box','off')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',16)
%     exportgraphics(gca,strcat('Figures/Polar_map_fi_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model, '.png'),'Resolution',300)
    exportgraphics(gca,strcat('Figures/Polar_map_fi_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model, '.jpeg'),'Resolution',300)

    %     saveas(gca,strcat('Figures/Polar_map_fi_y_idx_',num2str((yy)),'_fi_',num2str(fi)),'png')


    %% PDF on Polar map (eigen vectors)
    % Prepare histogram
    figure;
    polarhistogram(theta_map_1,n_bins_pdf, 'DisplayStyle','stairs','LineWidth',2, 'LineStyle','-','EdgeColor',[0 0.4470 0.7410]); %'Normalization', 'pdf'
    hold on
    polarhistogram(theta_map_2,n_bins_pdf,'DisplayStyle','stairs', 'LineWidth',2, 'LineStyle','-','EdgeColor',[0.4660 0.6740 0.1880]); %'Normalization', 'pdf'
    polarhistogram(theta_map_3,n_bins_pdf,'DisplayStyle','stairs', 'LineWidth',2, 'LineStyle','-','EdgeColor',[0.8500 0.3250 0.0980]); %'Normalization', 'pdf'
    %     thetalim([0 180])
    legend([{'$\theta_1$'},{'$\theta_2$'},{'$\theta_3$'}],'interpreter','latex','location','west','box','off')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',16)
%     saveas(gca,strcat('Figures/Polar_map_theta_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model),'png')
%     saveas(gca,strcat('Figures/Polar_map_theta_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model),'epsc')
%     exportgraphics(gca,strcat('Figures/Polar_map_theta_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model, '.png'),'Resolution',300)
    exportgraphics(gca,strcat('Figures/Polar_map_theta_y_idx_',num2str((yy)),'_fi_',num2str(fi), '_', SFS_model, '.jpeg'),'Resolution',300)



end