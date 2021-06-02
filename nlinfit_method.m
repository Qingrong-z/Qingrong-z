
%% mu_ref and T_ref
y1=ratio{1,1}(680:821);
x1=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta01=beta0_set{1};
[beta1,r,J,cov] = nlinfit(x1,y1,modelfun,beta01,options);
ci1 = nlparci(beta1,r,'covar',cov);
c21 = nlparci(beta1,r,'Jacobian',J);

y2=ratio{1,2}(680:821);
x2=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta02=beta0_set{2};
% beta02=[297;1.5;1.5];
[beta2,r,J,cov] = nlinfit(x2,y2,modelfun,beta02,options);
ci2 = nlparci(beta2,r,'covar',cov);
c22 = nlparci(beta2,r,'Jacobian',J);

y3=ratio{1,3}(680:821);
x3=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta03=beta0_set{3};
[beta3,r,J,cov] = nlinfit(x3,y3,modelfun,beta03,options);
ci3 = nlparci(beta3,r,'covar',cov);
c23 = nlparci(beta3,r,'Jacobian',J);

y4=ratio{1,4}(680:821);
x4=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta04=beta0_set{4};
[beta4,r,J,cov] = nlinfit(x4,y4,modelfun,beta04,options);
ci4 = nlparci(beta4,r,'covar',cov);
c24 = nlparci(beta4,r,'Jacobian',J);

y5=ratio{1,5}(680:821);
x5=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta05=beta0_set{5};
[beta5,r,J,cov] = nlinfit(x5,y5,modelfun,beta05,options);
ci5 = nlparci(beta5,r,'covar',cov);
c25 = nlparci(beta5,r,'Jacobian',J);

y6=ratio{1,6}(680:821);
x6=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta06=beta0_set{6};
[beta6,r,J,cov] = nlinfit(x6,y6,modelfun,beta06,options);
ci6 = nlparci(beta6,r,'covar',cov);
c26 = nlparci(beta6,r,'Jacobian',J);

y7=ratio{1,7}(680:821);
x7=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta07=beta0_set{7};
[beta7,r,J,cov] = nlinfit(x7,y7,modelfun,beta07,options);
ci7 = nlparci(beta7,r,'covar',cov);
c27 = nlparci(beta7,r,'Jacobian',J);

y8=ratio{1,8}(680:821);
x8=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta08=beta0_set{8};
[beta8,r,J,cov] = nlinfit(x8,y8,modelfun,beta08,options);
ci8 = nlparci(beta8,r,'covar',cov);
c28 = nlparci(beta8,r,'Jacobian',J);

y9=ratio{1,9}(680:821);
x9=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta09=beta0_set{9};
[beta9,r,J,cov] = nlinfit(x9,y9,modelfun,beta09,options);
ci9 = nlparci(beta9,r,'covar',cov);
c29 = nlparci(beta9,r,'Jacobian',J);

y10=ratio{1,10}(680:821);
x10=E(Einv(1.46):Einv(1.57));
modelfun=@(a,x) mu_T_ratio(a(1), a(2), a(3),E(Einv(1.46):Einv(1.57)),T_L,m_e,m_h,Eg,D);
options = statset('FunValCheck','off','DerivStep',10^-12,'robust','on');
beta010=beta0_set{10};
[beta10,r,J,cov] = nlinfit(x10,y10,modelfun,beta010,options);
ci10 = nlparci(beta10,r,'covar',cov);
c210 = nlparci(beta10,r,'Jacobian',J);

T_nfit1 = [beta1(1,1),beta2(1,1),beta3(1,1),beta4(1,1),beta5(1,1),beta6(1,1),beta7(1,1),beta8(1,1),beta9(1,1),beta10(1,1)];
mu_nfit1 = [beta1(2,1),beta2(2,1),beta3(2,1),beta4(2,1),beta5(2,1),beta6(2,1),beta7(2,1),beta8(2,1),beta9(2,1),beta10(2,1)];
mu_ref_nfit1 = [beta1(3,1),beta2(3,1),beta3(3,1),beta4(3,1),beta5(3,1),beta6(3,1),beta7(3,1),beta8(3,1),beta9(3,1),beta10(3,1)];


p_nift1 = polyfit(T_nfit1-T_L,P_abs,1); %Temperature fit
Q_nfit1 = p_nift1(1);
deltaT0_nift1 = p_nift1(2);
deltaT_nfit1data = deltaT0_nift1 + Q_nfit1*(T_nfit1-T_L);

p_nfit2 = polyfit(log(P_abs),mu_nfit1,1); %mu fit
mu_slope_nfit1 = p_nfit2(1);
mu_ref_nfit10 = p_nfit2(2);
mu_nfit1data = mu_ref_nfit10+mu_slope_nfit1*log(P_abs);
%%  Plots
%Mu_ref
figure
semilogx(P_abs, mu_ref_nfit1,'LineWidth',3)
ylabel('$\mu_{ref}$ (eV)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
% ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on


%Mu
figure
scatter(P_abs, mu_nfit1,200,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
hold on
plot(P_abs, mu_nfit1data,'--','color', colors(1,:),'LineWidth',2)
axis
ylabel('$\mu$ (eV)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([P_abs(1) P_abs(end)])
% ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'XScale','log')
set(gcf,'color','w')
box on


%T
figure
hold on
scatter(P_abs,T_nfit1-T_L,100,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
plot(deltaT_nfit1data,T_nfit1-T_L,'--','color', colors(2,:),'LineWidth', 2);
ylabel('$T-T_L$ (K)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on

%% Find the optimal mu_1

mu_ref_vec1=linspace(min(mu_ref_nfit1),max(mu_ref_nfit1));

for i2=1:length(mu_ref_vec1)
    for i1=2:length(Int)
        ratio_fit1{i2}=@(x) mu_T_ratio(x(1), x(2), mu_ref_vec1(i2),E(Einv(E_min):Einv(E_max)),T_L,m_e,m_h,Eg,D);
        ratio_opt1{i1,i2}=@(x) 1-abs(ratio_fit1{i2}(x)./ratio{i1}(Einv(E_min):Einv(E_max)));
        [x_sol1{i1,i2},res_sol1(i1,i2)]=lsqnonlin(ratio_opt1{i1,i2},x_init,x_min,x_max);
    end
end

res_sol_sum1=(sum(res_sol1,1));
[min_res1,idx_res1]=min(res_sol_sum1);

figure
for i1=2:length(Int)
    semilogy(mu_ref_vec1,res_sol1(i1,:),'linewidth',2)
    leg{i1-1}=['Ratio ' num2str(i1)];
    hold on
end
semilogy(mu_ref_vec1,res_sol_sum1,'color','black','linewidth',3)
leg{length(Int)}='Total';
xlabel('$\mu_{ref}$ (eV)','Interpreter','Latex')
ylabel('Residual error','Interpreter','Latex')
box on
xlim([min(mu_ref_vec1) max(mu_ref_vec1)])
ylim([0 inf])
set(gca,'Fontsize',12)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
legend(leg)
box on


mu_ref_bynfit=mu_ref_vec1(idx_res1);

for i1=2:length(Int)
    T_bynift(i1)=x_sol1{i1,idx_res1}(1);
    mu_bynift(i1)=x_sol1{i1,idx_res1}(2);
end
T_bynift(1)=T_L;
mu_bynift(1)=mu_ref_bynfit;