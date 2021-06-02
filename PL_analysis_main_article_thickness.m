clear all
close all

global E;

h=6.626e-34; %Planck's constant
h_ = h/(2*pi); %Reduced Planck's constant
c=2.998e8; %Speed of light
q=1.602e-19; % Elementary charge
k = 1.381e-23; %Boltzmann's constant
E_LO = 36.1e-3; %Energy of the LO phonon in GaAs (eV)
T_L = 296; % Room temperature


% % Get the value of the gap with the low intensity experience (no heat)
% Eg0 = 1.52; %eV at 0K
% alpha = 8.871*10^(-4); % eV/kelvin
% beta = 572; %Kelvin

%% Choose the sample

% day: '2019_07_05'
% sample: '1872_1_4' '1873_1_4' '1874_1_2'
% detector: 'vis'
% laser: 532

% day: '2019_09_16'
% sample: '1872_3_1' '1873_1_4' '1873_3_3' '1874_3_1'
% detector: 'vis'
% laser: 532

%day: '2019_11_28'
%sample: 'MBE2_3' 'MBE2_4' 'MBE2_6'
% detector: 'vis' 'IR'
% laser: 532, 638 (only for MBE2_3 & vis)

%day: '2019_12_10'
%sample: '1927_2_1' '1926_2_5' 'MBE1_2_1'
% detector: 'vis'
% laser: 532

%day: '2019_12_11'
%sample: '1872_1_4' '1873_1_4' '1874_1_2' '1926_2_5' '1927_2_1' 'MBE1_2_1'
% detector: 'vis'
% laser: 638

%day: '2020_01_15'
%sample: 'MBE1_3_2' 'MBE2_Q3_2_2'
%detector: 'vis' for MBE1_3_2, 'IR' for MBE2_Q3_2_2
%laser: 532, 638
day='2019_07_05';
sample='1873_1_4';
detector='vis';
laser=532;

[data, power, spot_surface, Int, E]=load_data_PL(day,sample,detector,laser);

%% Define parameters

switch sample
    case '1872_1_4' %GaAs 20nm Au mirror
        thickness=20e-7;
        barrier_thickness=85;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1872_wav_300_1000_habs_20_h_both_AlGaAs_84_Au_mirror.mat');
        A_GaAs_E_min1=A_GaAs_E;
        A_total_E_min1=A_total_E;
        Abs_layer_E_min1=Abs_layer_E;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1872_wav_300_1000_habs_20_h_both_AlGaAs_86_Au_mirror.mat');
        A_GaAs_E_plus1=A_GaAs_E;
        A_total_E_plus1=A_total_E;
        Abs_layer_E_plus1=Abs_layer_E;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1872_wav_300_1000_habs_20_h_both_AlGaAs_85_Au_mirror.mat');
        
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
        
    case '1874_1_2' %In0.1GaAs 20nm Au mirror
        thickness=20e-7;
        barrier_thickness=85;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1872_wav_300_1000_habs_20_h_both_AlGaAs_85_Au_mirror.mat');
    case '1872_3_1' %GaAs 20nm Ag mirror
        thickness=20e-7;
        barrier_thickness=91;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1872_wav_300_1000_habs_20_h_both_AlGaAs_91_Ag_mirror.mat');
    case '1873_3_3' %GaAs 200nm Ag mirror
        thickness=200e-7;
        barrier_thickness=92;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1873_wav_300_1000_habs_200_h_both_AlGaAs_92_Ag_mirror.mat');
    case '1874_3_1' %In0.1GaAs 20nm Ag mirror
        thickness=20e-7;
        barrier_thickness=91;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1872_wav_300_1000_habs_20_h_both_AlGaAs_91_Ag_mirror.mat');
    case '1926_2_5' %GaAs 50nm Au mirror
        thickness=50e-7;
        barrier_thickness=82;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1926_wav_300_1000_habs_50_h_both_AlGaAs_81_Au_mirror')
        A_GaAs_E_min1=A_GaAs_E;
        A_total_E_min1=A_total_E;
        Abs_layer_E_min1=Abs_layer_E;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1926_wav_300_1000_habs_50_h_both_AlGaAs_83_Au_mirror')
        A_GaAs_E_plus1=A_GaAs_E;
        A_total_E_plus1=A_total_E;
        Abs_layer_E_plus1=Abs_layer_E;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1926_wav_300_1000_habs_50_h_both_AlGaAs_82_Au_mirror')
        
    case '1927_2_1' %GaAs 100nm Au mirror
        thickness=100e-7;
        barrier_thickness=80;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1927_wav_300_1000_habs_100_h_both_AlGaAs_79_Au_mirror')
        A_GaAs_E_min1=A_GaAs_E;
        A_total_E_min1=A_total_E;
        Abs_layer_E_min1=Abs_layer_E;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1927_wav_300_1000_habs_100_h_both_AlGaAs_81_Au_mirror')
        A_GaAs_E_plus1=A_GaAs_E;
        A_total_E_plus1=A_total_E;
        Abs_layer_E_plus1=Abs_layer_E;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_1927_wav_300_1000_habs_100_h_both_AlGaAs_80_Au_mirror')
        
    case 'MBE1_2_1' %GaAs 20nm Au mirror (MBE)
        thickness = 20e-7; %(cm)
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_MBE1_wav_300_1000_habs_20_h_both_AlGaAs_72_Au_mirror')
    case 'MBE1_3_2'
        thickness = 20e-7; %(cm)
        barrier_thickness=88;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_MBE1_3_2_wav_300_1000_habs_20_h_both_AlGaAs_88_Au_mirror')
    case 'MBE2_Q3_2_2'
        thickness = 20e-7; %(cm)
        barrier_thickness=81;
        load('absorption\Abs_TMM_AlGaAs_GaAs_HCSC_ELO_MBE2_Q3_2_2_wav_300_1000_habs_20_h_both_AlGaAs_81_Au_mirror')
end


switch day
    case '2019_07_05'
        switch sample
            
            case '1872_1_4' %GaAs 20nm Au mirror
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[1 2 3 4 5];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.48;
                E_max = 1.50;
                E_ratio = 1.51;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
            case '1873_1_4' %GaAs 200nm Au mirror
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[1 2 3 12];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.53;
                E_max = 1.57;
                E_ratio = 1.51;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
            case '1874_1_2' %In0.1GaAs 20nm Au mirror
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[1 2 12];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.35;
                E_max = 1.38;
                E_ratio = 1.51;
                E_plot=[1.25 1.5];
                Eg = 1.31; %Band gap (eV)
        end
        
    case '2019_09_16'
        switch sample
            
            case '1873_1_4' %GaAs 200nm Au mirror
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[9 10];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.53;
                E_max = 1.57;
                E_ratio = 1.51;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
            case '1872_3_1' %GaAs 20nm Ag mirror
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[1 2 3 4];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.48;
                E_max = 1.5;
                E_ratio = 1.51;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
                
            case '1873_3_3' %GaAs 200nm Ag mirror
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[8 9 10];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.5;
                E_max = 1.6;
                E_ratio = 1.51;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
                
            case '1874_3_1' %InGaAs 20nm Ag mirror
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[7 8 9 10];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.35;
                E_max = 1.38;
                E_ratio = 1.51;
                E_plot=[1.25 1.5];
                Eg = 1.31; %Band gap (eV)
        end
    case '2019_11_28'
        switch sample
            case 'MBE2_3'
                thickness = 20e-7; %(cm)
                switch detector
                    case 'IR'
                        E_min = 1.17;
                        E_max = 1.2;
                        E_plot=[1 1.3];
                        Eg = 1.10; %Band gap (eV)
                    case 'vis'
                        E_min = 1.5;
                        E_max = 1.6;
                        E_plot=[1.35 1.6];
                        Eg = 1.424; %Band gap (eV)
                end
            case 'MBE2_4'
                thickness = 20e-7; %(cm)
                switch detector
                    case 'IR'
                        E_min = 1.17;
                        E_max = 1.2;
                        E_plot=[1 1.3];
                        Eg = 1.10; %Band gap (eV)
                    case 'vis'
                        E_min = 1.5;
                        E_max = 1.6;
                        E_plot=[1.35 1.6];
                        Eg = 1.424; %Band gap (eV)
                end
            case 'MBE2_6'
                thickness = 20e-7; %(cm)
                switch detector
                    case 'IR'
                        E_min = 1.17;
                        E_max = 1.2;
                        E_plot=[1 1.3];
                        Eg = 1.10; %Band gap (eV)
                    case 'vis'
                        E_min = 1.5;
                        E_max = 1.6;
                        E_plot=[1.35 1.6];
                        Eg = 1.424; %Band gap (eV)
                end
        end
    case '2019_12_10'
        switch sample
            case '1927_2_1'
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[7 8 9];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.53;
                E_max = 1.57;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
            case '1926_2_5'
                
                points_removed=[1 7];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.53;
                E_max = 1.57;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
            case 'MBE1_2_1'
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[6 7 8 9];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.53;
                E_max = 1.55;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
        end
        
    case '2019_12_11'
        switch sample
            case '1872_1_4'
                
                %removing points on which the model does not apply (lattice heating...)
                points_removed=[1 2 3 4];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.48;
                E_max = 1.50;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
            case '1926_2_5'
                
                points_removed=[];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.53;
                E_max = 1.57;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
            case '1927_2_1'
                
                points_removed=[];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.53;
                E_max = 1.57;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
            case '1873_1_4'
                
                points_removed=[1];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.53;
                E_max = 1.57;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
                
            case '1874_1_2'
                
                points_removed=[];
                Int(points_removed)=[];
                power(points_removed)=[]; %power is in W.cm^-2
                
                E_min = 1.35;
                E_max = 1.38;
                E_plot=[1.25 1.5];
                Eg = 1.31; %Band gap (eV)
                
            case 'MBE1_2_1'
                E_min = 1.49;
                E_max = 1.51;
                E_plot=[1.35 1.6];
                Eg = 1.424; %Band gap (eV)
        end
        
    case '2020_01_15'
        switch sample
            case 'MBE1_3_2'
                switch laser
                    case 532
                        %removing points on which the model does not apply (lattice heating...)
                        %                         points_removed=[7 8 9 10 11 12];
                        points_removed=[];
                        Int(points_removed)=[];
                        power(points_removed)=[]; %power is in W.cm^-2
                        
                        E_min = 1.50;
                        E_max = 1.60;
                        E_plot=[1.35 1.6];
                        Eg = 1.424; %Band gap (eV)
                        
                    case 638
                        %removing points on which the model does not apply (lattice heating...)
                        points_removed=[];
                        %                         points_removed=[];
                        Int(points_removed)=[];
                        power(points_removed)=[]; %power is in W.cm^-2
                        
                        E_min = 1.50;
                        E_max = 1.55;
                        E_plot=[1.35 1.6];
                        Eg = 1.424; %Band gap (eV)
                end
                
            case 'MBE2_Q3_2_2'
                switch laser
                    case 532
                        %removing points on which the model does not apply (lattice heating...)
                        points_removed=[];
                        Int(points_removed)=[];
                        power(points_removed)=[]; %power is in W.cm^-2
                        
                        E_min = 1.2;
                        E_max = 1.25;
                        E_plot=[1 1.3];
                        Eg = 1.10; %Band gap (eV)
                        
                    case 638
                        %removing points on which the model does not apply (lattice heating...)
                        points_removed=[1 2];
                        Int(points_removed)=[];
                        power(points_removed)=[]; %power is in W.cm^-2
                        
                        E_min = 1.2;
                        E_max = 1.25;
                        E_plot=[1 1.3];
                        Eg = 1.10; %Band gap (eV)
                end
        end
end

colors=lines(length(Int));
P_inc=power/(1000*spot_surface); % incident power (W.cm^-2)
E_laser=h*c/(q*laser*1e-9); %Energy of the laser (eV)

%% PL spectra

figure
for i=1:length(Int)
    semilogy(E,Int{i},'linewidth',2,'color',colors(i,:));
    hold on
end
xlim(E_plot)
% ylim([1e36 1e44])
xlabel('$E \: \mathrm{(eV)}$','Interpreter','Latex')
ylabel('$\phi_n$ (arb. units)','Interpreter','Latex')
box on
set(gca,'Fontsize',18)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');

%% Calculate the ratios and the temperature
T_raw(1) = T_L;
T_iter(1)=T_L;
delta_T_0_raw = 1;
delta_T_raw(1)=delta_T_0_raw;

for i1 = 1:length(Int)
    for i2 = i1:length(Int)
        
        % Basic ratios
        ratio{i1,i2}= real((Int{i2}./Int{i1}));
        fit_ratio = fit(E(Einv(E_min):Einv(E_max)), log(ratio{i1,i2}(Einv(E_min):Einv(E_max))),'poly1');
        coeffs_ratio = coeffvalues(fit_ratio);
        a1{i1,i2} = coeffs_ratio(1); %% unit = 1/eV
        a0{i1,i2} = coeffs_ratio(2); %% no unit mu2/(k*T2) - mu1/(k*T1)
        % Confidence intervals
        conf_int{i1,i2} = confint(fit_ratio,0.95);
        a1_sup{i1,i2} = conf_int{i1,i2}(2,1);
        delta_a1{i1,i2}=a1_sup{i1,i2}-a1{i1,i2};
        
        theo_log_ratio{i1,i2} = a1{i1,i2}*E+a0{i1,i2};
        %         ideality_factor(i1,i2) = (log(average_ratio{i1,i2}) - (1/T_raw(i1)-1/T_raw(i2))*q/k * mean(E(ref_ratio))) /log(power(i2)/power(i1));
        
        % Band filling terms
        log_BF{i1,i2} = real(log(ratio{i1,i2}) - theo_log_ratio{i1,i2});
        
        fit_ratio_min = fit(E(Einv(E_min-0.01):Einv(E_max-0.01)), log(ratio{i1,i2}(Einv(E_min-0.01):Einv(E_max-0.01))),'poly1');
        coeffs_ratio_min = coeffvalues(fit_ratio_min);
        a1_min{i1,i2} = coeffs_ratio_min(1); %% unit = 1/eV
        fit_ratio_plus = fit(E(Einv(E_min+0.01):Einv(E_max+0.01)), log(ratio{i1,i2}(Einv(E_min+0.01):Einv(E_max+0.01))),'poly1');
        coeffs_ratio_plus = coeffvalues(fit_ratio_plus);
        a1_plus{i1,i2} = coeffs_ratio_plus(1); %% unit = 1/eV
    end
end

% Temperature
for i=2:length(Int)
    T_raw(i) = 1/(1/T_raw(1)- k*a1{1,i}/q);
    T_plus_raw(i)=1/(1/T_raw(1)- k*a1_plus{1,i}/q);
    T_min_raw(i)=1/(1/T_raw(1)- k*a1_min{1,i}/q);
    delta_T_interval_raw(i)=max(abs(T_plus_raw(i)-T_raw(i)),abs(T_raw(i)-T_min_raw(i)));
    delta_T_slope_raw(i)=T_raw(i)^2*sqrt((delta_T_raw(1)/T_raw(1)^2)^2+(k/q)^2*delta_a1{1,i}^2);
    delta_T_raw(i)=delta_T_interval_raw(i)+delta_T_slope_raw(i);
end

%% Plot the PL ratios

% Compared to lowest intensity
figure
for i2=1:length(Int)
    semilogy(E, real((ratio{1,i2})),'color', colors(i2,:),'linewidth',2)
    hold on
    semilogy(E,exp(theo_log_ratio{1,i2}),'color', colors(i2,:))
end
xlim(E_plot)
xlabel('$E \: \mathrm{(eV)}$','Interpreter','Latex')
ylabel('$\phi_n/\phi_1$','Interpreter','Latex')
box on
set(gca,'Fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on

%% Absorbed power

A_barrier_E=Abs_layer_E(1,:)+Abs_layer_E(3,:); % absorption in the barriers
A_laser_GaAs=A_GaAs_E(round(interp1(E_A, (1:length(E_A)),E_laser))); % absorption in GaAs at laser wavelength
A_laser_GaAs_min1=A_GaAs_E_min1(round(interp1(E_A, (1:length(E_A)),E_laser))); % absorption in GaAs at laser wavelength
A_laser_GaAs_plus1=A_GaAs_E_plus1(round(interp1(E_A, (1:length(E_A)),E_laser))); % absorption in GaAs at laser wavelength

A_laser_barriers=A_barrier_E(round(interp1(E_A, (1:length(E_A)),E_laser))); % Absorption in barriers at laser wavelength
Barriers_to_GaAs=0.1; % Estimated fraction of carriers generated in the barriers that recombine in the absorber
A_laser=A_laser_GaAs+Barriers_to_GaAs*A_laser_barriers; % Equivalent absorptivity in the absorber
P_abs=A_laser*P_inc; %Equivalent absorbed power in the absorber
P_gen=(P_abs/thickness)*(1-Eg/E_laser);

% Errors
delta_power_rel=0.1; %relative error on power
delta_spot_surface_rel=0.2; %relative error on spot_surface
delta_A_GaAs=max(abs(A_laser_GaAs_plus1-A_laser_GaAs),abs(A_laser_GaAs-A_laser_GaAs_min1));
delta_A_laser_rel=(delta_A_GaAs+0.1*A_laser_barriers)/A_laser; %relative error on absorptivity
delta_P_abs_rel=sqrt(delta_power_rel^2+delta_spot_surface_rel^2+delta_A_laser_rel^2); %relative error on the absorbed power: systematic error
delta_P_gen_rel=delta_P_abs_rel;


%% Calculate Q_raw

[Q_fit_raw,delta_Q_fit_raw]=linfitxy([T_raw-T_L],[P_abs],[delta_T_raw],zeros(length(T_raw),1));
Q_raw=Q_fit_raw(1);
Q0=Q_fit_raw(2);
delta_Q_fit_raw_rel=delta_Q_fit_raw(1)/Q_fit_raw(1);
delta_Q_raw_rel=sqrt(delta_P_abs_rel^2+delta_Q_fit_raw_rel^2);
delta_Q_raw=Q_raw*delta_Q_raw_rel;



%% Modify Q_raw into Q

delta_T_0=sqrt(delta_T_raw(1)^2+(P_abs(1)/Q_raw)^2*(delta_P_abs_rel^2+delta_Q_raw_rel^2));
T(1)=T_L+P_abs(1)/Q_raw;
delta_T(1)=sqrt(delta_T_0^2+(P_abs(1)/Q_raw)^2*(delta_P_abs_rel^2+delta_Q_raw_rel^2));

for i=2:length(Int)
    T(i) = 1/(1/T(1)- k*a1{1,i}/q);
    T_plus(i)=1/(1/T(1)- k*a1_plus{1,i}/q);
    T_min(i)=1/(1/T(1)- k*a1_min{1,i}/q);
    delta_T_interval(i)=max(abs(T_plus(i)-T(i)),abs(T(i)-T_min(i)));
    delta_T_slope(i)=T(i)^2*sqrt((delta_T(1)/T(1)^2)^2+(k/q)^2*delta_a1{1,i}^2);
    delta_T(i)=delta_T_interval(i)+delta_T_slope(i);
end

% Q_fit=fit([0 T-T_L]',[0 P_abs]','poly1');
[Q_fit,delta_Q_fit]=linfitxy([0 T-T_L],[0 P_abs],[delta_T_0 delta_T],zeros(length(T_raw)+1,1));
Q=Q_fit(1);
delta_Q_fit_rel=delta_Q_fit(1)/Q_fit(1);
delta_Q_rel=sqrt(delta_P_abs_rel^2+delta_Q_fit_rel^2);
delta_Q=Q*delta_Q_rel;



% figure
% hold on
% errorbar(P_abs,T_raw-T_L,delta_T_raw,delta_T_raw,[],[], '+', 'markerSize',10,'color', colors(1,:),'LineWidth', 2);
% plot(Q_raw.*([0 T_raw]-T_L) + Q0,[0 T_raw]-T_L, 'color', colors(2,:),'LineWidth', 2);
% ylabel('$T-T_L$ (K)','Interpreter','Latex')
% xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
% box on
% xlim([0 P_abs(end)])
% % ylim([0 inf])
% set(gca,'Fontsize',12)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w')
% box on

figure
hold on
errorbar([0 P_abs],[0 T-T_L],[delta_T_0 delta_T],[delta_T_0 delta_T],[],[], '+', 'markerSize',10,'color', colors(1,:),'LineWidth', 2);
plot([0 Q.*(T-T_L)],[0 T-T_L], 'color', colors(2,:),'LineWidth', 2);
ylabel('$T-T_L$ (K)','Interpreter','Latex')
xlabel('$P_{abs} \: (\mathrm{W.cm^{-2})}$','Interpreter','Latex')
box on
xlim([0 P_abs(end)])
ylim([0 inf])
set(gca,'Fontsize',12)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w')
box on

%legend only
figure
for i=length(Int):-1:1
    semilogy(E(1),Int{i}(1),'linewidth',2,'color',colors(i,:));
    leg{i}=['P_{abs}=' num2str(round(P_abs(i))) ' W.cm^{-2}'];
    hold on
end
legend(flip(leg))
axis off
set(gca,'Fontsize',12)
set(gcf,'color','w');