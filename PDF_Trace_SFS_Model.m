function PDF_Trace_SFS_Model(x,y,z,idx_filt,Tau_xx_DNS_filt, Tau_xy_DNS_filt, Tau_xz_DNS_filt, Tau_yy_DNS_filt, Tau_yz_DNS_filt, Tau_zz_DNS_filt, ...
   Tau_kk_SFS,u_b,delta_h,fi)

%% Eigen Decomposition
% close all; clc

% Initialization
num_points_x = length(x(:,1,1));
num_points_y = length(x(1,:,1));
num_points_z = length(x(1,1,:));

% Y - positions for barycentric map
y_outer_vec     = [0.2 1 1.8]*delta_h;

% Cell array containing histogram data for each model and y-position
bar_tau_kk_norm = cell(cell(length(y_outer_vec),length(Tau_kk_SFS)+1));

for yy = 1:length(y_outer_vec)

    [~,idx_y]    = min(abs(y(1,:,1)-y_outer_vec(yy)));

    % Reset bar map
    bar_kk_map      = cell(cell(1,length(Tau_kk_SFS)+1));

    for nTrace = 1:length(Tau_kk_SFS)+1

 
            % Decomposition point by point
            for jj = idx_y:idx_y
                for ii = idx_filt:num_points_x-(idx_filt-1) % Only filtered points from y_outer
                    for kk = idx_filt:num_points_z-(idx_filt-1)


                        if nTrace <= 2 % SFS trace model

                            Tau_kk = Tau_kk_SFS{nTrace};

                            tau_filt_kk = Tau_kk(ii,jj,kk);

                        else % Filtered DNS

                            tau_ij_DNS_filt = [Tau_xx_DNS_filt(ii,jj,kk), Tau_xy_DNS_filt(ii,jj,kk), Tau_xz_DNS_filt(ii,jj,kk);
                                Tau_xy_DNS_filt(ii,jj,kk), Tau_yy_DNS_filt(ii,jj,kk), Tau_yz_DNS_filt(ii,jj,kk);
                                Tau_xz_DNS_filt(ii,jj,kk), Tau_yz_DNS_filt(ii,jj,kk), Tau_zz_DNS_filt(ii,jj,kk)];

                            tau_filt_kk = tau_ij_DNS_filt(1,1) + tau_ij_DNS_filt(2,2) + tau_ij_DNS_filt(3,3);

                        end

                        bar_kk_map{nTrace} = [bar_kk_map{nTrace}, tau_filt_kk];



                    end

                end
            end


    %% PDF magnitude (trace)
    bar_tau_kk_norm{yy,nTrace} = bar_kk_map{nTrace}./u_b^2;
    n_bins_pdf                 = 2*ceil(1 + log2(length(bar_kk_map{nTrace})));

    end

    %% PDF - Trace y/delta
    figure;
    histogram(bar_tau_kk_norm{yy,3},n_bins_pdf,'DisplayStyle','stairs', 'LineWidth',2, 'LineStyle','-','EdgeColor',[0 0.4470 0.7410]); %'Normalization', 'pdf'
    hold on
    histogram(bar_tau_kk_norm{yy,1},n_bins_pdf,'DisplayStyle','stairs', 'LineWidth',2, 'LineStyle','-','EdgeColor',[0.4660 0.6740 0.1880]); %'Normalization', 'pdf'
    histogram(bar_tau_kk_norm{yy,2},n_bins_pdf,'DisplayStyle','stairs', 'LineWidth',2, 'LineStyle','-','EdgeColor',[0.8500 0.3250 0.0980]); %'Normalization', 'pdf'
    xlabel('${\tau_{kk}/ {u_b}}^2$','interpreter','latex')
    ylabel('$PDF$','interpreter','latex')
    % legend([{'$y/\delta = 0.2$'},{'$y/\delta = 1.0$'},{'$y/\delta = 1.8$'}],'interpreter','latex','location','best')
    legend([{'$DNS filt$'},{'${SFS}_1$'},{'${SFS}_2$'}],'interpreter','latex','location','best')
    pbaspect([1 1 1])
    legend('Location','northeast','box','off')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
    if yy == 1
        xlim([-0.005 0.06])
    elseif yy == 2
        xlim([-0.001 0.02])
        xticklabels( {'0','0.005', '0.01', '0.015', '0.02'} )

    elseif yy == 3
        xlim([-0.005 0.06])
    end
%     saveas(gca,strcat('Figures/PDF_DNS_filtered_trace_y_idx_',num2str((yy)), '_fi_',num2str(fi), '_SFS_Comparison'),'png')
%     exportgraphics(gca,strcat('Figures/PDF_DNS_filtered_trace_y_idx_',num2str((yy)), '_fi_',num2str(fi), '_SFS_Comparison', '.png'),'Resolution',300)
%     exportgraphics(gca,strcat('Figures/PDF_DNS_filtered_trace_y_idx_',num2str((yy)), '_fi_',num2str(fi), '_SFS_Comparison', '.eps'),'Resolution',300)
    exportgraphics(gca,strcat('Figures/PDF_DNS_filtered_trace_y_idx_',num2str((yy)), '_fi_',num2str(fi), '_SFS_Comparison', '.jpeg'),'Resolution',300)



end



end