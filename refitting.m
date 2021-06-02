
% 
% x=E(Einv(1.46):Einv(1.57));
% modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
% options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
% beta0=[296;1.1;1.1];
% [beta,r,J,cov,mse] = nlinfit(x,ratio{1,2}(680:821),modelfun,beta0,options);
% ci = nlparci(beta,r,'covar',cov);
% c2 = nlparci(beta,r,'Jacobian',J);
% 
% %% mu_ref and T_ref
% y=ratio{1,9}(680:821);
% x=E(Einv(1.46):Einv(1.57));
% modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
% options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
% beta0=[296;1.17;1.148];
% beta = nlinfit(x,y,modelfun,beta0,options);
% ci = nlparci(beta,r,'covar',cov);
% c2 = nlparci(beta,r,'Jacobian',J);
%% Remove the second data set
% Int(2) = [];
% ratio(2) = [];
% T(2) = [];
% mu(2) = [];
% P_abs(2) = [];
% 
%% Feedback iteration
%     beta0=[295;1.2;1.1782];
% for i = 1:length(Int)
%     y = ratio{1,i}(Einv(1.46):Einv(1.57));
%     x = E(Einv(1.46):Einv(1.57));
%     modelfun = @(a,x) mu_T_ratio(a(1), a(2), a(3),x,T_L,m_e,m_h,Eg,D);
%     options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
%     [beta{i},r,J,cov,MSE(i)] = nlinfit(x,y,modelfun,beta0,options);
%     T_nfit(i) = beta{i}(1,1);
%     mu_nfit(i) = beta{i}(2,1);
%     mu_ref_nfit(i) = beta{i}(3,1);
%     ci{i} = nlparci(beta{i},r,'covar',cov);
%     cj{i} = nlparci(beta{i},r,'Jacobian',J);
%     beta0 = beta{i};
% end
% 
% p1 = polyfit(T_nfit-T_L,P_abs,1); %Temperature fit
% Q_nfit = p1(1);
% deltaT0_nift = p1(2);
% deltaT_nfitdata = deltaT0_nift + Q_nfit*(T_nfit-T_L);
% 
% p2 = polyfit(log(P_abs),mu_nfit,1); %mu fit
% mu_slope_nfit = p2(1);
% mu_ref_nfit0 = p2(2);
% mu_nfitdata = mu_ref_nfit0+mu_slope_nfit*log(P_abs);
% 
% figure
% title('iteration fit')
% hold on
% scatter(P_abs,T_nfit-T_L,100,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
% plot(deltaT_nfitdata,T_nfit-T_L,'--','color', colors(2,:),'LineWidth', 2);
% ylabel('$T-T_L$ (K)','Interpreter','Latex')
% xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
% box on
% xlim([0 P_abs(end)])
% ylim([0 inf])
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w')
% box on
% 
% figure
% title('iteration fit')
% semilogx(P_abs, mu_ref_nfit,'LineWidth',3)
% ylabel('$\mu_{ref}$ (eV)','Interpreter','Latex')
% xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
% box on
% xlim([0 P_abs(end)])
% % ylim([0 inf])
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w')
% box on
% 
% figure
% title('iteration fit')
% scatter(P_abs, mu_nfit,200,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
% hold on
% plot(P_abs, mu_nfitdata,'--','color', colors(1,:),'LineWidth',2)
% axis
% ylabel('$\mu$ (eV)','Interpreter','Latex')
% xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
% box on
% xlim([P_abs(1) P_abs(end)])
% % ylim([0 inf])
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca,'XScale','log')
% set(gcf,'color','w')
% box on

%% Initial set by previous data
% for i = 1:length(Int)
%     y = ratio{1,i}(Einv(1.46):Einv(1.57));
%     x = E(Einv(1.46):Einv(1.57));
%     modelfun = @(a,x) mu_T_ratio(a(1), a(2), a(3),x,T_L,m_e,m_h,Eg,D);
%     options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
%     beta0_set{i} = [T(i);mu(i);mu_ref];
%     [beta_set{i},r,J,cov] = nlinfit(x,y,modelfun,beta0_set{i},options);
%     T_nfit_set(i) = beta_set{i}(1,1);
%     mu_nfit_set(i) = beta_set{i}(2,1);
%     mu_ref_nfit_set(i) = beta_set{i}(3,1);
%     ci_set{i} = nlparci(beta_set{i},r,'covar',cov);
%     cj_set{i} = nlparci(beta_set{i},r,'Jacobian',J);
% end
% 
% p_set1 = polyfit(T_nfit_set-T_L,P_abs,1); %Temperature fit
% Q_set_nfit = p_set1(1);
% deltaT0_set_nift = p_set1(2);
% deltaT_set_nfitdata = deltaT0_set_nift + Q_set_nfit*(T_nfit_set-T_L);
% 
% p_set2 = polyfit(log(P_abs),mu_nfit_set,1); %mu fit
% mu_slope_set_nfit = p_set2(1);
% mu_ref_set_nfit0 = p_set2(2);
% mu_set_nfitdata = mu_ref_set_nfit0+mu_slope_set_nfit*log(P_abs);
% % 
% figure
% title('initial set')
% hold on
% scatter(P_abs,T_nfit_set-T_L,100,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
% plot(deltaT_set_nfitdata,T_nfit_set-T_L,'--','color', colors(2,:),'LineWidth', 2);
% ylabel('$T-T_L$ (K)','Interpreter','Latex')
% xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
% box on
% xlim([0 P_abs(end)])
% ylim([0 inf])
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w')
% box on
% 
% figure
% title('initial set')
% semilogx(P_abs, mu_ref_nfit_set,'LineWidth',3)
% ylabel('$\mu_{ref}$ (eV)','Interpreter','Latex')
% xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
% box on
% xlim([0 P_abs(end)])
% % ylim([0 inf])
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w')
% box on
% 
% figure
% title('initial set')
% scatter(P_abs, mu_nfit_set,200,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
% hold on
% plot(P_abs, mu_set_nfitdata,'--','color', colors(1,:),'LineWidth',2)
% axis
% ylabel('$\mu$ (eV)','Interpreter','Latex')
% xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
% box on
% xlim([P_abs(1) P_abs(end)])
% % ylim([0 inf])
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca,'XScale','log')
% set(gcf,'color','w')
% box on

%% Replace mu as delta_mu + mu_ref
% % Initialize to find the range of mu_1 values
% ratio_dmu_fit_init=@(a) mu_T_ratio(a(1), a(2), a(3),E(Einv(E_min):Einv(E_max)),T_L,m_e,m_h,Eg,D);
% for i1 = 1:length(Int)
%     ratio{i1}= real(Int{i1}./Int{1});
%     ratio_opt_init_dmu{i1}=@(a,E) abs(1-ratio_dmu_fit_init(a)./ratio{i1}(Einv(E_min):Einv(E_max)));
%     
%     y = ratio_opt_init_dmu{i1}(Einv(E_min):Einv(E_max));
%     x = E(Einv(E_min):Einv(E_max));
% %     modelfun = @(a,x) dmu_T_ratio(a(1), a(2), a(3),x,T_L,m_e,m_h,Eg,D);
%     options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
%     beta0_dmu_init{i1} = [T_L;0.3;0.9];
%     [beta_dmu_init{i1},r,J,cov] = nlinfit(x,y,ratio_opt_init_dmu{i1},beta0_dmu_init{i1},options);
%     
%     T_nfit_dmu_init(i1)=beta_dmu_init{i1}(1);
%     dmu_nfit_dmu_init(i1)=beta_dmu_init{i1}(2);
%     mu_ref_nfit_dmu_init(i1)=beta_dmu_init{i1}(3);
% end
% 
% T_nfit_dmu_init(1)=T_L;
% mu_nfit_dmu_init(1)=NaN;
% mu_ref_nfit_dmu_init(1)=NaN;


% Find optimal mu_ref
% mu_ref_vec_dmu=linspace(min(mu_ref_nfit_dmu_init),max(mu_ref_nfit_dmu_init));
% 
% for i2=1:length(mu_ref_vec_dmu)
%     for i1=2:length(Int)
%         ratio_fit_dmu{i2}=@(a,x) dmu_T_ratio(a(1), a(2), mu_ref_vec_dmu(i2),E(Einv(E_min):Einv(E_max)),T_L,m_e,m_h,Eg,D);
%         ratio_opt_dmu{i1,i2}=@(x) 1-abs(ratio_fit{i2}(x)./ratio{i1}(Einv(E_min):Einv(E_max)));
%         [x_sol{i1,i2},res_sol(i1,i2)]=lsqnonlin(ratio_opt{i1,i2},x_init,x_min,x_max);
%     end
% end
% 
% res_sol_sum=(sum(res_sol,1));
% [min_res,idx_res]=min(res_sol_sum);
% 


% Initialize to find the range of mu_1 values

for i = 2:length(Int)
    y = ratio{1,i}(Einv(E_min):Einv(E_max));
    x = E(Einv(E_min):Einv(E_max));
    modelfun = @(a,x) dmu_T_ratio(a(1), a(2), a(3),x,T_L,m_e,m_h,Eg,D);
    options = statset('FunValCheck','off','DerivStep',10^-13,'robust','on');
    beta0_dmu_init{i} = [T(i);mu(i)-mu_ref;mu_ref];
%     beta0_dmu_init{i} = [T(i);mu(i)-mu_ref_init(i);mu_ref_init(i)];
    [beta_dmu_init{i},r,J,cov] = nlinfit(x,y,modelfun,beta0_dmu_init{i},options);
    T_nfit_dmu_init(i) = beta_dmu_init{i}(1,1);
    delta_mu_init(i) = beta_dmu_init{i}(2,1);
    mu_nfit_dmu_init(i) = beta_dmu_init{i}(2,1)+beta_dmu_init{i}(3,1);
    mu_ref_nfit_dmu_init(i) = beta_dmu_init{i}(3,1);
    ci_dmu_init{i} = nlparci(beta_dmu_init{i},r,'covar',cov);
    cj_dmu_init{i} = nlparci(beta_dmu_init{i},r,'Jacobian',J);
    mu_ref_inv_var(i) = 4/(ci_dmu_init{i}(3,2)-ci_dmu_init{i}(3,1));
    mu_top(i) =  mu_ref_nfit_dmu_init(i).*(mu_ref_inv_var(i).^2);
end
 mu_ref_nfit_dmu_init(1)=NaN;
 mu_ref_nift_dmu = sum(mu_top)./sum((mu_ref_inv_var.^2));
 var_mu_ref_nfit_dmu = 2.*sqrt(1./sum((mu_ref_inv_var.^2)));


figure
title('f(deltamu);initial set')
semilogx(P_abs, mu_ref_nfit_dmu_init,'LineWidth',3)
ylabel('$\mu_{ref}$ (eV)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
% ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on
 
T_nfit_dmu_init(1)=T_L;
%  delta_mu_init(1)=NaN;
%  mu_nfit_dmu_init(1)=NaN;
%  mu_ref_nfit_dmu_init(1)=NaN;

% Refit with obtained optimal mu_ref
for i = 1:length(Int)
    y = ratio{1,i}(Einv(E_min):Einv(E_max));
    x = E(Einv(E_min):Einv(E_max));
    modelfun = @(a,x) dmu_T_ratio(a(1), a(2), mu_ref_nift_dmu,x,T_L,m_e,m_h,Eg,D);
    options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
    beta0_dmu{i} = [T_nfit_dmu_init(i);delta_mu_init(i)];
    [beta_dmu{i},r,J,cov] = nlinfit(x,y,modelfun,beta0_dmu{i},options);
    T_nfit_dmu(i) = beta_dmu{i}(1,1);
    delta_mu(i) = beta_dmu{i}(2,1);
    mu_nfit_dmu(i) = beta_dmu{i}(2,1)+mu_ref_nift_dmu;
    ci_dmu{i} = nlparci(beta_dmu{i},r,'covar',cov);
    cj_dmu{i} = nlparci(beta_dmu{i},r,'Jacobian',J);
    T_error_dmu(i) = (ci_dmu{i}(1,2)-ci_dmu{i}(1,1))/2;
    delta_mu_error(i) = (ci_dmu{i}(2,2)-ci_dmu{i}(2,1))/2;
end

% Linear fit
p_dmu1 = polyfit(T_nfit_dmu-T_L,P_abs,1); %Temperature fit
Q_dmu_nfit = p_dmu1(1);
deltaT0_dmu_nift = p_dmu1(2);
deltaT_dmu_nfitdata = deltaT0_dmu_nift + Q_dmu_nfit*(T_nfit_dmu-T_L);

p_dmu2 = polyfit(log(P_abs),delta_mu,1); %mu fit
mu_slope_dmu_nfit = p_dmu2(1);
mu_ref_dmu_nfit0 = p_dmu2(2);
mu_dmu_nfitdata = mu_ref_dmu_nfit0+mu_slope_dmu_nfit*log(P_abs);

% Plots
figure
hold on
errorbar(P_abs,T_nfit_dmu-T_L,10*T_error_dmu,10*T_error_dmu,[],[], '+', 'markerSize',10,'color', colors(1,:),'LineWidth', 3);
plot(deltaT_dmu_nfitdata,T_nfit_dmu-T_L,'--','color', colors(2,:),'LineWidth', 2);
yyaxis left
ylabel('$T-T_L$ (K)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on
yyaxis right
ylabel('$Error$ (K)','Interpreter','Latex')
ylim([0 5])



figure
errorbar(P_abs,delta_mu,50*delta_mu_error,50*delta_mu_error,[],[], '+', 'markerSize',10,'color', colors(1,:),'LineWidth', 3);
hold on
plot(P_abs, mu_dmu_nfitdata,'--','color', colors(2,:),'LineWidth',2)
axis
yyaxis left
ylabel('$\Delta\mu$ (eV)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([P_abs(1) P_abs(end)])
ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'XScale','log')
set(gcf,'color','w')
box on
yyaxis right
ylabel('$Error$ (eV)','Interpreter','Latex')
ylim([0 0.004])
%% Log ratio fitting
% Initialize to find the range of mu_1 values
for i = 2:length(Int)
    y = log(ratio{1,i}(Einv(E_min):Einv(E_max)));
    x = E(Einv(E_min):Einv(E_max));
    modelfun = @(a,x) log(dmu_T_ratio(a(1), a(2), a(3),x,T_L,m_e,m_h,Eg,D));
    options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
    beta0_log_set_init{i} = [T(i);mu(i)-mu_ref;mu_ref];
%     beta0_log_set_init{i} = [T(i);mu(i);mu_ref_init(i)];
    [beta_log_set_init{i},r,J,cov] = nlinfit(x,y,modelfun,beta0_log_set_init{i},options);
    T_nfit_log_set_init(i) = beta_log_set_init{i}(1,1);
    dmu_nfit_log_set_init(i) = beta_log_set_init{i}(2,1);
    mu_ref_nfit_log_set_init(i) = beta_log_set_init{i}(3,1);
    ci_log_set_init{i} = nlparci(beta_log_set_init{i},r,'covar',cov);
    cj_log_set_init{i} = nlparci(beta_log_set_init{i},r,'Jacobian',J);
    mu_ref_inv_var_log(i) = 4/(ci_log_set_init{i}(3,2)-ci_log_set_init{i}(3,1));
    mu_top_log(i) =   mu_ref_nfit_log_set_init(i).*(mu_ref_inv_var_log(i).^2);
end
mu_ref_nfit_log_set_init(1)=NaN;
mu_ref_nift_log = sum(mu_top_log)./sum((mu_ref_inv_var_log.^2));
var_mu_ref_nfit_log = 2.*sqrt(1./sum((mu_ref_inv_var_log.^2)));



figure
title('log ratio')
semilogx(P_abs, mu_ref_nfit_log_set_init,'LineWidth',3)
ylabel('$\mu_{ref}$ (eV)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
% ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on

T_nfit_log_set_init(1)=T_L;
% dmu_nfit_log_set_init(1)=0;
% Refit with obtained optimal mu_ref
for i = 1:length(Int)
    y = log(ratio{1,i}(Einv(E_min):Einv(E_max)));
    x = E(Einv(E_min):Einv(E_max));
    modelfun = @(a,x) log(dmu_T_ratio(a(1), a(2), mu_ref_nift_log,x,T_L,m_e,m_h,Eg,D));
    options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
    beta0_log_set{i} = [T_nfit_log_set_init(i);dmu_nfit_log_set_init(i)];
    [beta_log_set{i},r,J,cov] = nlinfit(x,y,modelfun,beta0_log_set{i},options);
    T_nfit_log_set(i) = beta_log_set{i}(1,1);
    dmu_nfit_log_set(i) = beta_log_set{i}(2,1);
    ci_log_set{i} = nlparci(beta_log_set{i},r,'covar',cov);
    cj_log_set{i} = nlparci(beta_log_set{i},r,'Jacobian',J);
    T_error_log(i) = (ci_log_set{i}(1,2)-ci_log_set{i}(1,1))/2;
    dmu_error_log(i) = (ci_log_set{i}(2,2)-ci_log_set{i}(2,1))/2;
end

p_log_set1 = polyfit(T_nfit_log_set-T_L,P_abs,1); %Temperature fit
Q_log_set_nfit = p_log_set1(1);
deltaT0_log_set_nift = p_log_set1(2);
deltaT_log_set_nfitdata = deltaT0_log_set_nift + Q_log_set_nfit*(T_nfit_log_set-T_L);

p_log_set2 = polyfit(log(P_abs),dmu_nfit_log_set,1); %mu fit
dmu_slope_log_set_nfit = p_log_set2(1);
dmu_ref_log_set_nfit0 = p_log_set2(2);
dmu_log_set_nfitdata = dmu_ref_log_set_nfit0+dmu_slope_log_set_nfit*log(P_abs);

figure
errorbar(P_abs,T_nfit_log_set-T_L,10*T_error_log,10*T_error_log,[],[], '+', 'markerSize',10,'color', colors(1,:),'LineWidth', 3);
% scatter(P_abs, T_nfit_log_set-T_L,200,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
hold on
plot(deltaT_log_set_nfitdata,T_nfit_log_set-T_L,'--','color', colors(2,:),'LineWidth', 2);
yyaxis left
ylabel('$T-T_L$ (K)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on
yyaxis right
ylabel('$Error$ (K)','Interpreter','Latex')
ylim([0 5])


figure
errorbar(P_abs,dmu_nfit_log_set,62.5*dmu_error_log,62.5*dmu_error_log,[],[], '+', 'markerSize',10,'color', colors(1,:),'LineWidth', 3);
% scatter(P_abs, dmu_nfit_log_set,200,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
hold on
plot(P_abs, dmu_log_set_nfitdata,'--','color', colors(2,:),'LineWidth',2)
axis
yyaxis left
ylabel('$\Delta\mu$ (eV)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([P_abs(1) P_abs(end)])
% ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'XScale','log')
set(gcf,'color','w')
box on
yyaxis right
yyaxis right
ylabel('$Error$ (eV)','Interpreter','Latex')
ylim([0 0.004])
%% Plots
figure
hold on
plot(P_abs, mu_dmu_nfitdata,'--','color', colors(1,:),'LineWidth',2)
errorbar(P_abs,delta_mu,62.5*delta_mu_error,62.5*delta_mu_error,[],[], '+', 'markerSize',10,'color', colors(1,:),'LineWidth', 3);
hold on
plot(P_abs, dmu_log_set_nfitdata,'--','color', colors(2,:),'LineWidth',2)
errorbar(P_abs,dmu_nfit_log_set,62.5*dmu_error_log,62.5*dmu_error_log,[],[], '+', 'markerSize',10,'color', colors(2,:),'LineWidth', 3);
axis
ylabel('$\mu$ (eV)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([P_abs(1) P_abs(end)])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'XScale','log')
set(gcf,'color','w')
legend('linear \Delta\mu','linear \Delta\mu error','log\Delta\mu','log\Delta\mu error')
box on

figure
hold on
plot(deltaT_dmu_nfitdata,T_nfit_dmu-T_L,'--','color', colors(1,:),'LineWidth', 2);
errorbar(P_abs,T_nfit_dmu-T_L,10*T_error_dmu,10*T_error_dmu,[],[], '+', 'markerSize',10,'color', colors(1,:),'LineWidth', 3);
hold on
plot(deltaT_log_set_nfitdata,T_nfit_log_set-T_L,'--','color', colors(2,:),'LineWidth', 2);
errorbar(P_abs,T_nfit_log_set-T_L,10*T_error_log,10*T_error_log,[],[], '+', 'markerSize',10,'color', colors(2,:),'LineWidth', 3);
axis
ylabel('$T-T_L$ (K)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
legend('linear \Delta\mu','linear \Delta\mu error','log\Delta\mu','log\Delta\mu error')
box on


