clearvars
close all

global E;

h=6.626e-34; %Planck's constant
h_ = h/(2*pi); %Reduced Planck's constant
c=2.998e8; %Speed of light
q=1.602e-19; % Elementary charge
k = 1.381e-23; %Boltzmann's constant
E_LO = 36.1e-3; %Energy of the LO phonon in GaAs (eV)
T_L = 296; % Room temperature

m_e=0.067; %Effective mass of electrons, relative to the free electron mass
m_h=0.47; %Effective mass of holes, relative to the free electron mass
D=3; %Dimensionality of bulk


% % Get the value of the gap with the low intensity experience (no heat)
% Eg0 = 1.52; %eV at 0K
% alpha = 8.871*10^(-4); % eV/kelvin
% beta = 572; %Kelvin

%% Choose the sample

day='2019_07_05';
sample='1873_1_4';
detector='vis';
laser=532;

[data, power, spot_surface, Int, E]=load_data_PL(day,sample,detector,laser);

%% Define parameters

switch sample
        
    case '1873_1_4' %GaAs 200nm Au mirror
        thickness=200e-7;
        barrier_thickness=93;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1873_wav_300_1000_habs_200_h_both_AlGaAs_92_Au_mirror.mat');
        A_GaAs_E_min1=A_GaAs_E;
        A_total_E_min1=A_total_E;
        Abs_layer_E_min1=Abs_layer_E;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1873_wav_300_1000_habs_200_h_both_AlGaAs_94_Au_mirror.mat');
        A_GaAs_E_plus1=A_GaAs_E;
        A_total_E_plus1=A_total_E;
        Abs_layer_E_plus1=Abs_layer_E;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1873_wav_300_1000_habs_200_h_both_AlGaAs_93_Au_mirror.mat');
        %find the optimal, uncertainty of the absorpiton
         
end

%
switch day
    case '2019_07_05'
        switch sample
                           
            case '1873_1_4' %GaAs 200nm Au mirror
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.46;
                E_max = 1.57;
                E_ratio = 1.51;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)

        end
end

colors=lines(length(Int));
P_inc=power/(1000*spot_surface); % incident power (W.cm^-2)
E_laser=h*c/(q*laser*1e-9); %Energy of the laser (eV)

%% PL spectra





%% Calculate the ratios with band filling
%Parameter intialization and limits
% T=zeros(size(Int));
% mu=zeros(size(Int));
% mu_ref=zeros(size(Int));

x_init=[T_L, 1.2, 1.2];
x_min=[T_L-10, 0.8, 0.8];
x_max=[500, Eg, Eg];

%% Initialize to find the range of mu_1 values

ratio_fit_init=@(x) mu_T_ratio(x(1), x(2), x(3),E(Einv(E_min):Einv(E_max)),T_L,m_e,m_h,Eg,D);

for i1 = 1:length(Int)
    ratio{i1}= real(Int{i1}./Int{1});
    ratio_opt_init{i1}=@(x) abs(1-ratio_fit_init(x)./ratio{i1}(Einv(E_min):Einv(E_max)));
    [x_sol_init{i1},res_sol_init(i1)]=lsqnonlin(ratio_opt_init{i1},x_init,x_min,x_max);
    T_init(i1)=x_sol_init{i1}(1);
    mu_init(i1)=x_sol_init{i1}(2);
    mu_ref_init(i1)=x_sol_init{i1}(3);
end

T_init(1)=T_L;
mu_init(1)=NaN;
mu_ref_init(1)=NaN;


% %Plot the PL ratios
% figure
% for i1=1:length(Int)
%     semilogy(E, ratio{i1},'--','color', colors(i1,:),'linewidth',2)
%     hold on
%     semilogy(E,mu_T_ratio(T_init(i1), mu_init(i1), mu_ref_init(i1),E,T_L,m_e,m_h,Eg,D),'color', colors(i1,:))
%     hold on
% end
% xlim(E_plot)
% xlabel('$E \: \mathrm{(eV)}$','Interpreter','Latex')
% ylabel('$\phi_n/\phi_1$','Interpreter','Latex')
% box on
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w')
% box on


%% Find the optimal mu_1

mu_ref_vec=linspace(min(mu_ref_init),max(mu_ref_init));

for i2=1:length(mu_ref_vec)
    for i1=2:length(Int)
        ratio_fit{i2}=@(x) mu_T_ratio(x(1), x(2), mu_ref_vec(i2),E(Einv(E_min):Einv(E_max)),T_L,m_e,m_h,Eg,D);
        ratio_opt{i1,i2}=@(x) 1-abs(ratio_fit{i2}(x)./ratio{i1}(Einv(E_min):Einv(E_max)));
        [x_sol{i1,i2},res_sol(i1,i2)]=lsqnonlin(ratio_opt{i1,i2},x_init,x_min,x_max);
    end
end

res_sol_sum=(sum(res_sol,1));
[min_res,idx_res]=min(res_sol_sum);




mu_ref=mu_ref_vec(idx_res);

for i1=2:length(Int)
    T(i1)=x_sol{i1,idx_res}(1);
    mu(i1)=x_sol{i1,idx_res}(2);
end
T(1)=T_L;
mu(1)=mu_ref;

%% Plot the PL ratios



%% Absorbed power

%% Calculate T and mu slopes

switch sample
    case '1873_1_4'
        T_interv=[3:length(T)];
        mu_interv=[1:5];
    
end

T_fit=polyfit(T(T_interv)-T_L,P_abs(T_interv),1);
Q=T_fit(1);
deltaT0=T_fit(2);

deltaTfitdata=deltaT0+Q*(T-T_L);

mu_fit=polyfit(log(P_abs(mu_interv)),mu(mu_interv),1);
mu_slope=mu_fit(1);
mu_0=mu_fit(2);
mufitdata=mu_0+mu_slope*log(P_abs);
ideality_coeff=mu_slope*q/(k*T_L);

%% Plots

figure
semilogx(P_abs, mu_ref_init,'LineWidth',3)
ylabel('$\mu_{ref}$ (eV)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
% ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on


figure
scatter(P_abs, mu,200,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
hold on
plot(P_abs, mufitdata,'--','color', colors(1,:),'LineWidth',2)
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

figure
hold on
scatter(P_abs,T-T_L,100,'x','MarkerEdgeColor', colors(1,:),'LineWidth',3);
plot(deltaTfitdata,T-T_L,'--','color', colors(2,:),'LineWidth', 2);
ylabel('$T-T_L$ (K)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
ylim([0 inf])
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on

% figure
% semilogx(P_abs, mu-mu_ref,'LineWidth',3)
% % hold on
% % semilogx(P_abs,[0 mu_init(2:end)-mu_ref_init(2:end)],'LineWidth',3)
% ylabel('$\mu$ increase (eV)','Interpreter','Latex')
% xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
% box on
% xlim([0 P_abs(end)])
% ylim([0 inf])
% set(gca,'Fontsize',16)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w')
% box on


%% Legend only
% figure
% for i=length(Int):-1:1
%     semilogy(E(1),Int{i}(1),'linewidth',2,'color',colors(i,:));
%     leg{i}=['P_{abs}=' num2str(round(P_abs(i))) ' W.cm^{-2}'];
%     hold on
% end
% legend(flip(leg))
% axis off
% set(gca,'Fontsize',12)
% set(gcf,'color','w');
