%% HIGH PRESSURE - RATIO DENSITY SCALE OVER KOLMOGOROV SCALE
% Load results
DATA        = readtable('Data/output_data_1D_HighPressure_200x1x1_2Pc_test.csv');

% Settings from HPCFS
bSolver           = 'Real';                    % Define Ideal or Real Gas Model
% Substance
Substance         = 'N2';                      % Fluid selected substance
HP_model          = 'HighPressure';                % Transport coefficients model: 'Constant', 'LowPressure', 'HighPressure'
% Fluid properties
[~, T_c, P_c, ~, ~, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Substance);
% PengRobinson
[a,b,R,dadT,d2adT2,NASA_coefficients] = PengRobinson(T, Substance);

% Index only unique Y and Z
Index_Z = unique(DATA.Z); Index_Z = find(DATA.Z == Index_Z(2));
Index = Index_Z;
x = unique(DATA.X);
y = unique(DATA.Y);
z = unique(DATA.Z);
[X,Y,Z] = meshgrid(x,y,z);

[dx,~,~] = CentralDerivative_d1_2ndOrder(X);
[~,dy,~] = CentralDerivative_d1_2ndOrder(Y);
[~,~,dz] = CentralDerivative_d1_2ndOrder(Z);

fi     = 2;

% Variables
P     = reshape(DATA.P,[length(y),length(x),length(z)]);
T     = reshape(DATA.T,[length(y),length(x),length(z)]);
rho   = reshape(DATA.rho,[length(y),length(x),length(z)]);

% Filter variable
P_filt    = FilterFields_1D_CDLF_Box(P,fi,dx);
P_filt_avg = P_filt(2,:,2);


rho_filt  = FilterFields_1D_CDLF_Box(rho,fi,dx);
T_filt    = FilterFields_1D_CDLF_Box(T,fi,dx);

P_LES     = Calculate_P_PengRobinson(rho_filt,T_filt,Substance);
P_LES_avg = P_LES(2,:,2);


alpha_6_1D  = P_filt - P_LES;


%% Taylor expansion model
J_t2.dcv_dT = 0;
J_t         = Jacobian_thermodynamics(bSolver, rho_filt,T_filt,P_filt, J_t2,Substance);
d2P_drhodT  = J_t.d2P_drhodT;

rhoT_filt    = FilterFields_1D_CDLF_Box(rho_filt.*T_filt,fi,dx);
delta_P      = 0.5*d2P_drhodT.*(rho_filt.*T_filt - rhoT_filt);
P_Taylor     = P_LES + delta_P;
P_Taylor_avg = P_Taylor(2,:,2);

%% ILA model
P_LES_filt      = FilterFields_1D_CDLF_Box(P_LES,fi,dx);

rho_filt_filt   = FilterFields_1D_CDLF_Box(rho_filt,fi,dx);
T_filt_filt     = FilterFields_1D_CDLF_Box(T_filt,fi,dx);

P_filt_LES      = Calculate_P_PengRobinson(rho_filt_filt,T_filt_filt,Substance);

T_filt_FG               = FilterFields_1D_CDLF_Box(T_filt,2*fi,dx);
T_filt_FG_filt          = FilterFields_1D_CDLF_Box(T_filt_FG,fi,dx);
T_filt_FG_filt_FG       = FilterFields_1D_CDLF_Box(T_filt_FG_filt,2*fi,dx);
rho_filt_FG             = FilterFields_1D_CDLF_Box(rho_filt,2*fi,dx);
rho_filt_FG_filt        = FilterFields_1D_CDLF_Box(rho_filt_FG,fi,dx);
rho_filt_FG_filt_FG     = FilterFields_1D_CDLF_Box(rho_filt_FG_filt,2*fi,dx);

P_filt_FG_LES           = Calculate_P_PengRobinson(rho_filt_FG,T_filt_FG,Substance);
P_filt_FG_LES_filt      = FilterFields_1D_CDLF_Box(P_filt_FG_LES,fi,dx);
P_filt_FG_LES_filt_FG   = FilterFields_1D_CDLF_Box(P_filt_FG_LES_filt,2*fi,dx);

P_filt_FG_filt_FG_LES   = Calculate_P_PengRobinson(rho_filt_FG_filt_FG,T_filt_FG_filt_FG,Substance);


P_LES_filt_FG           = FilterFields_1D_CDLF_Box(P_LES_filt,2*fi,dx);
P_filt_LES_FG           = FilterFields_1D_CDLF_Box(P_filt_LES,2*fi,dx);

M_P = (P_filt_FG_LES_filt_FG - P_filt_FG_filt_FG_LES) - (P_LES_filt_FG - P_filt_LES_FG);


P_LES_FG   = FilterFields_1D_CDLF_Box(P_LES,2*fi,dx);

L_P         = P_LES_FG - P_filt_FG_LES;

% Coefficient C_p 
num_points_x = length(P(:,2,2));
num_points_y = length(P(2,:,2));
num_points_z = length(P(2,2,:));

C_P = L_P.*M_P./(M_P.*M_P);
C_P(C_P<0) = 0;


P_ILA = P_LES + C_P.*(P_LES_filt - P_filt_LES);
P_ILA_avg = P_ILA(2,:,2);


idx = 3*fi + 1; delta_h = 1;
% EOS
figure; hold on; box on
% set(gca,'yscale','log')
semilogy(x(idx:end-idx-1)/delta_h,(P_LES_avg(idx:end-idx-1))/P_c,'linewidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
hold on
semilogy(x(idx:end-idx-1)/delta_h,(P_filt_avg(idx:end-idx-1))/P_c,'linewidth',2, 'LineStyle','-','color',[0.4660 0.6740 0.1880])
semilogy(x(idx:end-idx-1)/delta_h,(P_Taylor_avg(idx:end-idx-1))/P_c,'linewidth',2, 'LineStyle','-.','color',[0.3010 0.7450 0.9330])
semilogy(x(idx:end-idx-1)/delta_h,(P_ILA_avg(idx:end-idx-1))/P_c,'linewidth',2, 'LineStyle','--','color',[0.8500 0.3250 0.0980])

% plot(y(2,idx:end-idx-1,2)/delta_h,alpha_6_avg{1}(idx:end-idx-1),'linewidth',2)
xlabel('${y/\delta}$','interpreter','latex')
ylabel('${P/P_c}$','interpreter','latex')
legend([{'$P(\overline{\rho},\breve{T})$'},{'$\overline{P(\rho,T)}$'},{'${P(\overline{\rho},\breve{T})} + {\delta_P}$'},{'${P_{ILA}}$'}],'interpreter','latex','location','best')
pbaspect([1.8 1 1])
legend('Location','southeast','box','off','NumColumns', 3)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 1]); ylim([1.999 2.001])
% saveas(gca,strcat('Figures/Unclosed_EOS_', num2str(fi), 'xDelta_Comparison_1D_test'),'png')
exportgraphics(gca,strcat('Figures/Unclosed_EOS_',num2str(fi), 'xDelta_Comparison_1D_test','.png'),'Resolution',300)
exportgraphics(gca,strcat('Figures/Unclosed_EOS_',num2str(fi), 'xDelta_Comparison_1D_test','.jpeg'),'Resolution',300)
