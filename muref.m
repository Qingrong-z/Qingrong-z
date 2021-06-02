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
%% Calculate the ratios with band filling
x_init=[T_L, 1.2, 1.2];
x_min=[T_L-10, 0.8, 0.8];
x_max=[500, Eg, Eg];
%% Initialize mu_ref
ratio_fit_init=@(x,E) mu_T_ratio(x(1), x(2), x(3),E(Einv(E_min):Einv(E_max)),T_L,m_e,m_h,Eg,D);

for i1 = 1:length(Int)
    ratio{i1}= real(Int{i1}./Int{1});
    beta0 = [300;1;1];
    beta{i1} = nlinfit(E,ratio{i1},mu_T_ratio(T, mu, mu_ref,E,T_L,m_e,m_h,Eg,D),beta0);
end
