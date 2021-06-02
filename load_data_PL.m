function [data, power, spot_surface, Int_cal, E]=load_data_PL(day,sample,detector,laser)

h=6.626e-34; %Planck's constant
c=2.998e8; %Speed of light
q=1.602e-19; % Elementary charge
delimiter = '\t';

load('Amaury_general\filters.mat'); %transmission ratio of the filters
direction_name = ['..\Amaury Setup\' day '\data\'];

switch day
    case '2019_07_05'
        ratio_power=0.497; %ratio between real power measured at sample position and where we actually measure it
        spot_surface = 1.15e-5; % surface of the laser spot (cm^2)
    case '2019_09_16'
        ratio_power=0.54; %ratio between real power measured at sample position and where we actually measure it
        spot_surface = 1.15e-5; % surface of the laser spot (cm^2)
    case '2019_11_28'
        switch laser
            case 532
                ratio_power=0.56; %ratio between real power measured at sample position and where we actually measure it
            case 638
                ratio_power=0.05; %Not measured
        end
        spot_surface=1.15e-5; % surface of the laser spot (cm^2)
        
    case '2019_12_10'
        ratio_power=0.56; %ratio between real power measured at sample position and where we actually measure it
        spot_surface = 1.15e-5; % surface of the laser spot (cm^2)
    case '2019_12_11'
        ratio_power=0.105;
        spot_surface=1.15e-5;
        
    case '2020_01_15'
        switch laser
            case 532
                ratio_power=0.53;
            case 638
                ratio_power=0.12;
        end
        spot_surface=1.15e-5;
    case '2020_06_18'
        ratio_power=0.52;
        spot_surface=1.15e-5;
end

switch detector
    case 'vis'
        formatSpec = '%f%f%f%[^\n\r]';
    case 'IR'
        formatSpec = '%f%f%[^\n\r]';
        startRow = 18;
        endRow = 529;
end

%% Generating the calibration vector 

filename = 'Amaury_general\089430034_HL-3 plus-CAL_cc_20190617_VIS(e).lmp';
fileID = fopen(filename,'r');%Open file for reading
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose('all');

% Takes the values
lambda_lamp = dataArray{:,1};
Int_lamp = dataArray{:,2}; %in 礧.cm^-2.nm^-1
E_lamp=flipud(h*c./(q*lambda_lamp*1e-9));%元素上下翻转

% Transforms the units to allow absolute calibration
Int_lamp=flipud(Int_lamp*1e7); %% data in W.m^-2.m^-1
Int_lamp = Int_lamp *h*c./(q*E_lamp).^2; %% data in s^-1.m^-2
Int_lamp = Int_lamp./(q*E_lamp); %% data in s^-1.m^-2.J^-1
Int_lamp = Int_lamp/(2*pi); %% data in s^-1.m^-2.J^-1.sr^-1 (2 pi because emission in a half sphere)

switch day
    case '2019_07_05'
        calib_file{1}.name='ocean_optics_50ms';
        calib_file{1}.acq=100;
        calib_file{1}.time=0.1;
        calib_file{1}.filter=1;
        
        calib_file{2}.name='thorlabs_fiber';
        calib_file{2}.acq=100;
        calib_file{2}.time=0.5;
        calib_file{2}.filter=1;
        
        calib_file{3}.name='thorlabs_sample';
        calib_file{3}.acq=5;
        calib_file{3}.time=20;
        calib_file{3}.filter=1;
        
    case '2019_09_16'
        calib_file{1}.name='oceanoptics_20190914';
        calib_file{1}.acq=100;
        calib_file{1}.time=0.1;
        calib_file{1}.filter=1;
        
        calib_file{2}.name='thorlabs_fiber';
        calib_file{2}.acq=50;
        calib_file{2}.time=0.1;
        calib_file{2}.filter=1;
        
        calib_file{3}.name='thorlabs_sample';
        calib_file{3}.acq=50;
        calib_file{3}.time=0.1;
        calib_file{3}.filter=1;
        
    case '2019_11_28'
        calib_file{1}.name='oceanoptics_20190914';
        calib_file{1}.acq=100;
        calib_file{1}.time=0.1;
        calib_file{1}.filter=1;
        
        calib_file{2}.name='thorlabs_fiber_NE20';
        calib_file{2}.acq=50;
        calib_file{2}.time=0.2;
        calib_file{2}.filter=NE20_exp_trans;
        
        calib_file{3}.name='thorlabs_sample_NE30';
        calib_file{3}.acq=50;
        calib_file{3}.time=1;
        calib_file{3}.filter=NE30_exp_trans;
        
    case '2019_12_10'
        
        calib_file{1}.name='oceanoptics_20190914';
        calib_file{1}.acq=100;
        calib_file{1}.time=0.1;
        calib_file{1}.filter=1;
        
        calib_file{2}.name='thorlabs_fiber';
        calib_file{2}.acq=50;
        calib_file{2}.time=0.15;
        calib_file{2}.filter=NE20_exp_trans;
        
        calib_file{3}.name='thorlabs_sample';
        calib_file{3}.acq=50;
        calib_file{3}.time=0.5;
        calib_file{3}.filter=NE20_exp_trans;
        
    case '2019_12_11'
        
        calib_file{1}.name='ocean_optics';
        calib_file{1}.acq=50;
        calib_file{1}.time=1;
        calib_file{1}.filter=NE20_exp_trans;
        
        calib_file{2}.name='thorlabs_fiber';
        calib_file{2}.acq=50;
        calib_file{2}.time=0.3;
        calib_file{2}.filter=NE20_exp_trans;
        
        calib_file{3}.name='thorlabs_sample';
        calib_file{3}.acq=50;
        calib_file{3}.time=1;
        calib_file{3}.filter=NE20_exp_trans;
        
    case '2020_01_15'
        
        calib_file{1}.name='ocean_optics';
        calib_file{1}.acq=50;
        calib_file{1}.time=1;
        calib_file{1}.filter=NE20_exp_trans;
        
        switch detector
            case 'vis'
                calib_file{2}.name='thorlabs_fiber';
                calib_file{2}.acq=50;
                calib_file{2}.time=0.15;
                calib_file{2}.filter=NE20_exp_trans;
                
                calib_file{3}.name='thorlabs_sample';
                calib_file{3}.acq=50;
                calib_file{3}.time=0.3;
                calib_file{3}.filter=NE20_exp_trans;
        
            case 'IR'
                calib_file{2}.name='thorlabs_fiber_nir';
                calib_file{2}.acq=10;
                calib_file{2}.time=2;
                calib_file{2}.filter=1;
                
                calib_file{3}.name='thorlabs_sample_nir';
                calib_file{3}.acq=10;
                calib_file{3}.time=5;
                calib_file{3}.filter=1;
        end
        
    case '2020_06_18'
        
        calib_file{1}.name='ocean_optics';
        calib_file{1}.acq=50;
        calib_file{1}.time=1;
        calib_file{1}.filter=NE20_exp_trans;
        
        calib_file{2}.name='thorlabs_fiber';
        calib_file{2}.acq=50;
        calib_file{2}.time=0.15;
        calib_file{2}.filter=NE20_exp_trans;
        
        calib_file{3}.name='thorlabs_sample';
        calib_file{3}.acq=50;
        calib_file{3}.time=0.3;
        calib_file{3}.filter=NE20_exp_trans;
   
end

for i=1:length(calib_file)
    
    for condition = ["bright", "dark"]
        %Finds the file
        filename = strcat(direction_name, calib_file{1}.name);
        if strcmp(condition, 'dark') == 1
            filename = strcat(filename,'_dark');
        end
        
        %Processes the file
        filename = strcat(filename,'.txt');
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
        fclose('all');
        lambda = dataArray{:,1};
        
        %Which object to fill with the data
        if strcmp(condition, 'dark') == 1
            Int_dark{i} = zeros(length(lambda),1);
            for j=2:size(dataArray,2)-1
                Int_dark{i} = Int_dark{i}(:)+dataArray{:,j};
            end
        else
            Int{i} = zeros(length(lambda),1);
            for j=2:size(dataArray,2)-1
                Int{i} = Int{i}+dataArray{:,j};
            end
        end
    end
    
    E=flipud(h*c./(q*lambda*1e-9));
    Int_abs{i}=flipud(Int{i}-Int_dark{i})./(calib_file{i}.acq*calib_file{i}.time.*calib_file{i}.filter);
    Int_abs{i}=filloutliers(abs(Int_abs{i}),'linear','movmedian',5,'ThresholdFactor',0.5);
    
end

ocean_optics_abs = Int_abs{1};
thorlabs_fiber_abs = Int_abs{2};
thorlabs_sample_abs = Int_abs{3};

Int_lamp_real = interp1(E_lamp,Int_lamp,E, 'spline');

ratio_real_ocean_lamp = Int_lamp_real./ocean_optics_abs; %ocean optics real / ocean optics measured when the lamp is before the fiber
ratio_optics = thorlabs_fiber_abs./thorlabs_sample_abs; %calculates the losses in the mirror between the sample and the 2nd fiber (with thorlabs lamp)
calib = ratio_real_ocean_lamp .* ratio_optics;


%% loading the sample into 'data'
switch day
    case '2019_07_05'
        switch sample
            case '1872_1_4' %GaAs 20nm Au mirror
                
                PM1_to_PM2 = 1.05; % multiplication factor PM1*PM1_av=PM2
                PM1_av=(1/2)*(1+PM1_to_PM2); %the average power,  when PM1 is measured
                PM2_av=(1/2)*(1+1/PM1_to_PM2); %the average power,  when PM1 is measured
                
                data{1}.pow = 0.372*ratio_power*PM1_av ;
                data{1}.acq = 40;
                data{1}.time = 1.5;
                data{1}.filter = 1;
                data{1}.name = 'a';
                
                data{2}.pow = 1.22*ratio_power*PM1_av ;
                data{2}.acq = 40;
                data{2}.time = 0.4;
                data{2}.filter = 1;
                data{2}.name = 'b';
                
                data{3}.pow = 3.67*ratio_power*PM1_av ;
                data{3}.acq = 40;
                data{3}.time = 1;
                data{3}.filter = NE10_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 13.9*ratio_power*PM2_av ;
                data{4}.acq = 40;
                data{4}.time = 0.7;
                data{4}.filter = NE20_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =26.8*ratio_power*PM2_av ;
                data{5}.acq = 40;
                data{5}.time = 1;
                data{5}.filter = NE30_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 52.5*ratio_power*PM2_av ;
                data{6}.acq = 40;
                data{6}.time = 0.5;
                data{6}.filter = NE30_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 101.5*ratio_power*PM2_av ;
                data{7}.acq = 40;
                data{7}.time = 1;
                data{7}.filter = NE40_exp_trans;
                data{7}.name = 'g';
                
                data{8}.pow = 201*ratio_power*PM2_av ;
                data{8}.acq = 40;
                data{8}.time = 0.4;
                data{8}.filter = NE40_exp_trans;
                data{8}.name = 'h';
                
                data{9}.pow = 284*ratio_power*PM2_av ;
                data{9}.acq = 40;
                data{9}.time = 0.2;
                data{9}.filter = NE40_exp_trans;
                data{9}.name = 'i';
                
                data{10}.pow = 400*ratio_power*PM2_av ;
                data{10}.acq = 40;
                data{10}.time = 1.2;
                data{10}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{10}.name = 'j';
                
                data{11}.pow = 561*ratio_power*PM2_av ;
                data{11}.acq = 40;
                data{11}.time = 0.9;
                data{11}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{11}.name = 'k';
                
                data{12}.pow = 798*ratio_power*PM2_av ;
                data{12}.acq = 40;
                data{12}.time = 0.5;
                data{12}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{12}.name = 'l';
                
            case '1873_1_4' %GaAs 200nm Au mirror
                
                PM1_to_PM2 = 1.05; % multiplication factor PM1*PM1_av=PM2
                PM1_av=(1/2)*(1+PM1_to_PM2); %the average power,  when PM1 is measured
                PM2_av=(1/2)*(1+1/PM1_to_PM2); %the average power,  when PM1 is measured
                
                data{1}.pow = 0.332*ratio_power*PM1_av ;
                data{1}.acq = 40;
                data{1}.time = 1.8;
                data{1}.filter = NE20_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 1.078*ratio_power*PM1_av ;
                data{2}.acq = 40;
                data{2}.time = 0.4;
                data{2}.filter = NE20_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 3.3*ratio_power*PM1_av ;
                data{3}.acq = 40;
                data{3}.time = 0.4;
                data{3}.filter = NE30_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 10.5*ratio_power*PM1_av ;
                data{4}.acq = 40;
                data{4}.time = 0.4;
                data{4}.filter = NE40_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =27.1*ratio_power*PM1_av ;
                data{5}.acq = 40;
                data{5}.time = 2.8;
                data{5}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 51.3*ratio_power*PM1_av ;
                data{6}.acq = 40;
                data{6}.time = 1.4;
                data{6}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 99*ratio_power*PM1_av ;
                data{7}.acq = 40;
                data{7}.time = 0.6;
                data{7}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{7}.name = 'g';
                
                data{8}.pow = 202*ratio_power*PM2_av ;
                data{8}.acq = 40;
                data{8}.time = 0.3;
                data{8}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{8}.name = 'h';
                
                data{9}.pow = 285*ratio_power*PM2_av ;
                data{9}.acq = 40;
                data{9}.time = 1.7;
                data{9}.filter = NE40_exp_trans.*NE20_exp_trans.*NE10_exp_trans;
                data{9}.name = 'i';
                
                data{10}.pow = 400*ratio_power*PM2_av ;
                data{10}.acq = 40;
                data{10}.time = 1.2;
                data{10}.filter = NE40_exp_trans.*NE20_exp_trans.*NE10_exp_trans;
                data{10}.name = 'j';
                
                data{11}.pow = 563*ratio_power*PM2_av ;
                data{11}.acq = 40;
                data{11}.time = 1;
                data{11}.filter = NE40_exp_trans.*NE20_exp_trans.*NE10_exp_trans;
                data{11}.name = 'k';
                
                data{12}.pow = 801*ratio_power*PM2_av ;
                data{12}.acq = 40;
                data{12}.time = 0.7;
                data{12}.filter = NE40_exp_trans.*NE20_exp_trans.*NE10_exp_trans;
                data{12}.name = 'l';
                
            case '1874_1_2' %InGaAs 20nm Au mirror
                
                PM1_to_PM2 = 1.05; % multiplication factor PM1*PM1_av=PM2
                PM1_av=(1/2)*(1+PM1_to_PM2); %the average power,  when PM1 is measured
                PM2_av=(1/2)*(1+1/PM1_to_PM2); %the average power,  when PM1 is measured
                
                data{1}.pow = 0.34*ratio_power*PM1_av ;
                data{1}.acq = 40;
                data{1}.time = 1;
                data{1}.filter = 1;
                data{1}.name = 'a';
                
                data{2}.pow = 1.1*ratio_power*PM1_av ;
                data{2}.acq = 40;
                data{2}.time = 0.3;
                data{2}.filter = 1;
                data{2}.name = 'b';
                
                data{3}.pow = 3.34*ratio_power*PM1_av ;
                data{3}.acq = 40;
                data{3}.time = 1.5;
                data{3}.filter = NE20_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 10.7*ratio_power*PM1_av ;
                data{4}.acq = 40;
                data{4}.time = 0.4;
                data{4}.filter = NE20_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =27.1*ratio_power*PM1_av ;
                data{5}.acq = 40;
                data{5}.time = 0.6;
                data{5}.filter = NE30_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 51.3*ratio_power*PM1_av ;
                data{6}.acq = 40;
                data{6}.time = 0.2;
                data{6}.filter = NE30_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 98.6*ratio_power*PM1_av ;
                data{7}.acq = 40;
                data{7}.time = 0.5;
                data{7}.filter = NE40_exp_trans;
                data{7}.name = 'g';
                
                data{8}.pow = 201*ratio_power*PM2_av ;
                data{8}.acq = 40;
                data{8}.time = 2.4;
                data{8}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{8}.name = 'h';
                
                data{9}.pow = 285*ratio_power*PM2_av ;
                data{9}.acq = 40;
                data{9}.time = 1.5;
                data{9}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{9}.name = 'i';
                
                data{10}.pow = 400*ratio_power*PM2_av ;
                data{10}.acq = 40;
                data{10}.time = 1;
                data{10}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{10}.name = 'j';
                
                data{11}.pow = 563*ratio_power*PM2_av ;
                data{11}.acq = 40;
                data{11}.time = 0.8;
                data{11}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{11}.name = 'k';
                
                data{12}.pow = 800*ratio_power*PM2_av ;
                data{12}.acq = 40;
                data{12}.time = 0.5;
                data{12}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{12}.name = 'l';
                
        end
        
    case '2019_09_16'
        switch sample
            case '1873_1_4' %GaAs 200nm Au mirror
                
                PM1_to_PM2 = 1.17; % multiplication factor PM1*PM1_av=PM2
                PM1_av=(1/2)*(1+PM1_to_PM2); %the average power,  when PM1 is measured
                PM2_av=(1/2)*(1+1/PM1_to_PM2); %the average power,  when PM1 is measured
                
                data{1}.pow = 1.12*ratio_power*PM1_av ;
                data{1}.acq = 50;
                data{1}.time = 0.3;
                data{1}.filter = NE20_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 3.36*ratio_power*PM1_av ;
                data{2}.acq = 50;
                data{2}.time = 0.3;
                data{2}.filter = NE30_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 11.2*ratio_power*PM1_av ;
                data{3}.acq = 50;
                data{3}.time = 0.5;
                data{3}.filter = NE30_exp_trans.*NE10_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 28.4*ratio_power*PM1_av ;
                data{4}.acq = 50;
                data{4}.time = 0.7;
                data{4}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =52.2*ratio_power*PM1_av ;
                data{5}.acq = 50;
                data{5}.time = 0.3;
                data{5}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 112*ratio_power*PM1_av ;
                data{6}.acq = 50;
                data{6}.time = 0.5;
                data{6}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 220*ratio_power*PM1_av ;
                data{7}.acq = 50;
                data{7}.time = 0.3;
                data{7}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{7}.name = 'g';
                
                data{8}.pow = 449*ratio_power*PM2_av ;
                data{8}.acq = 50;
                data{8}.time = 0.5;
                data{8}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{8}.name = 'h';
                
                data{9}.pow = 601*ratio_power*PM2_av ;
                data{9}.acq = 50;
                data{9}.time = 0.5;
                data{9}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{9}.name = 'i';
                
                data{10}.pow = 750*ratio_power*PM2_av ;
                data{10}.acq = 50;
                data{10}.time = 0.5;
                data{10}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{10}.name = 'j';
                
            case '1872_3_1'
                
                PM1_to_PM2 = 1.14;
                PM1_av=(1/2)*(1+PM1_to_PM2); %the average power,  when PM1 is measured
                PM2_av=(1/2)*(1+1/PM1_to_PM2); %the average power,  when PM1 is measured
                
                data{1}.pow = 1.13*ratio_power*PM1_av ;
                data{1}.acq = 50;
                data{1}.time = 2;
                data{1}.filter = 1;
                data{1}.name = 'a';
                
                data{2}.pow = 3.43*ratio_power*PM1_av ;
                data{2}.acq = 50;
                data{2}.time = 0.5;
                data{2}.filter = 1;
                data{2}.name = 'b';
                
                data{3}.pow = 11.2*ratio_power*PM1_av ;
                data{3}.acq = 50;
                data{3}.time = 2;
                data{3}.filter = NE10_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 24.1*ratio_power*PM1_av ;
                data{4}.acq = 50;
                data{4}.time = 0.5;
                data{4}.filter = NE10_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =54.4*ratio_power*PM1_av ;
                data{5}.acq = 50;
                data{5}.time = 0.2;
                data{5}.filter = NE10_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 110*ratio_power*PM1_av ;
                data{6}.acq = 50;
                data{6}.time = 0.2;
                data{6}.filter = NE20_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 222*ratio_power*PM1_av ;
                data{7}.acq = 50;
                data{7}.time = 0.2;
                data{7}.filter = NE30_exp_trans;
                data{7}.name = 'g';
                
                data{8}.pow = 450*ratio_power*PM2_av ;
                data{8}.acq = 50;
                data{8}.time = 0.5;
                data{8}.filter = NE40_exp_trans;
                data{8}.name = 'h';
                
                data{9}.pow = 601*ratio_power*PM2_av ;
                data{9}.acq = 50;
                data{9}.time = 0.2;
                data{9}.filter = NE40_exp_trans;
                data{9}.name = 'i';
                
                data{10}.pow = 748*ratio_power*PM2_av ;
                data{10}.acq = 50;
                data{10}.time = 0.2;
                data{10}.filter = NE40_exp_trans;
                data{10}.name = 'j';
                
                data{11}.pow = 900*ratio_power*PM2_av ;
                data{11}.acq = 50;
                data{11}.time = 0.2;
                data{11}.filter = NE40_exp_trans;
                data{11}.name = 'k';
                
            case '1873_3_3'
                
                PM1_to_PM2 = 1.12;
                PM1_av=(1/2)*(1+PM1_to_PM2); %the average power,  when PM1 is measured
                PM2_av=(1/2)*(1+1/PM1_to_PM2); %the average power,  when PM1 is measured
                
                data{1}.pow = 1.19*ratio_power*PM1_av ;
                data{1}.acq = 50;
                data{1}.time = 0.3;
                data{1}.filter = NE20_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 3.54*ratio_power*PM1_av ;
                data{2}.acq = 50;
                data{2}.time = 0.3;
                data{2}.filter = NE30_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 11.7*ratio_power*PM1_av ;
                data{3}.acq = 50;
                data{3}.time = 0.3;
                data{3}.filter = NE40_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 24.8*ratio_power*PM1_av ;
                data{4}.acq = 50;
                data{4}.time = 0.7;
                data{4}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =54.4*ratio_power*PM1_av ;
                data{5}.acq = 50;
                data{5}.time = 0.3;
                data{5}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 112*ratio_power*PM1_av ;
                data{6}.acq = 50;
                data{6}.time = 0.3;
                data{6}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 221*ratio_power*PM1_av ;
                data{7}.acq = 50;
                data{7}.time = 0.7;
                data{7}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{7}.name = 'g';
                
                data{8}.pow = 450*ratio_power*PM2_av ;
                data{8}.acq = 50;
                data{8}.time = 0.7;
                data{8}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{8}.name = 'h';
                
                data{9}.pow = 600*ratio_power*PM2_av ;
                data{9}.acq = 50;
                data{9}.time = 0.7;
                data{9}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{9}.name = 'i';
                
                data{10}.pow = 750*ratio_power*PM2_av ;
                data{10}.acq = 50;
                data{10}.time = 0.7;
                data{10}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{10}.name = 'j';
                
            case '1874_3_1'
                
                PM1_to_PM2 = 1.14;
                PM1_av=(1/2)*(1+PM1_to_PM2); %the average power,  when PM1 is measured
                PM2_av=(1/2)*(1+1/PM1_to_PM2); %the average power,  when PM1 is measured
                
                data{1}.pow = 1.16*ratio_power*PM1_av ;
                data{1}.acq = 50;
                data{1}.time = 2;
                data{1}.filter = NE20_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 3.51*ratio_power*PM1_av ;
                data{2}.acq = 50;
                data{2}.time = 0.5;
                data{2}.filter = NE20_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 11.5*ratio_power*PM1_av ;
                data{3}.acq = 50;
                data{3}.time = 0.5;
                data{3}.filter = NE30_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 24.6*ratio_power*PM1_av ;
                data{4}.acq = 50;
                data{4}.time = 1;
                data{4}.filter = NE40_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =53.9*ratio_power*PM1_av ;
                data{5}.acq = 50;
                data{5}.time = 0.5;
                data{5}.filter = NE40_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 113*ratio_power*PM1_av ;
                data{6}.acq = 50;
                data{6}.time = 0.15;
                data{6}.filter = NE40_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 221*ratio_power*PM1_av ;
                data{7}.acq = 50;
                data{7}.time = 0.5;
                data{7}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{7}.name = 'g';
                
                data{8}.pow = 450*ratio_power*PM2_av ;
                data{8}.acq = 50;
                data{8}.time = 0.15;
                data{8}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{8}.name = 'h';
                
                data{9}.pow = 600*ratio_power*PM2_av ;
                data{9}.acq = 50;
                data{9}.time = 0.15;
                data{9}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{9}.name = 'i';
                
                data{10}.pow = 750*ratio_power*PM2_av ;
                data{10}.acq = 50;
                data{10}.time = 0.25;
                data{10}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{10}.name = 'j';
        end
    case '2019_11_28'
        switch sample
            case 'MBE2_3'
                switch detector
                    case 'vis'
                        switch laser
                            case 638
                                data{1}.pow = 7.2*ratio_power ;
                                data{1}.acq = 50;
                                data{1}.time = 3;
                                data{1}.filter = NE30_exp_trans;
                                data{1}.name = 'a';
                                
                            case 532
                                data{1}.pow = 7.17*ratio_power ;
                                data{1}.acq = 50;
                                data{1}.time = 0.5;
                                data{1}.filter = NE30_exp_trans;
                                data{1}.name = 'a';
                                
                                data{2}.pow = 20.0*ratio_power ;
                                data{2}.acq = 50;
                                data{2}.time = 0.5;
                                data{2}.filter = NE40_exp_trans;
                                data{2}.name = 'b';
                                
                                data{3}.pow = 68.7*ratio_power ;
                                data{3}.acq = 50;
                                data{3}.time = 0.15;
                                data{3}.filter = NE40_exp_trans;
                                data{3}.name = 'c';
                                
                                data{4}.pow = 200*ratio_power;
                                data{4}.acq = 50;
                                data{4}.time = 0.5;
                                data{4}.filter = NE40_exp_trans.*NE20_exp_trans;
                                data{4}.name = 'd';
                                
                                data{5}.pow =450*ratio_power;
                                data{5}.acq = 50;
                                data{5}.time = 0.5;
                                data{5}.filter = NE40_exp_trans.*NE20_exp_trans;
                                data{5}.name = 'e';
                                
                                data{6}.pow = 700*ratio_power;
                                data{6}.acq = 50;
                                data{6}.time = 0.5;
                                data{6}.filter = NE40_exp_trans.*NE20_exp_trans;
                                data{6}.name = 'f';
                        end
                        
                        
                    case 'IR'
                        data{1}.pow = 7.17*ratio_power ;
                        data{1}.acq = 5;
                        data{1}.time = 30;
                        data{1}.filter = 1;
                        data{1}.name = 'a';
                        
                        data{2}.pow = 20.0*ratio_power ;
                        data{2}.acq = 5;
                        data{2}.time = 10;
                        data{2}.filter = 1;
                        data{2}.name = 'b';
                        
                        data{3}.pow = 69.7*ratio_power ;
                        data{3}.acq = 5;
                        data{3}.time = 3;
                        data{3}.filter = 1;
                        data{3}.name = 'c';
                        
                        data{4}.pow = 200*ratio_power;
                        data{4}.acq = 5;
                        data{4}.time = 1;
                        data{4}.filter = 1;
                        data{4}.name = 'd';
                        
                        data{5}.pow =450*ratio_power;
                        data{5}.acq = 5;
                        data{5}.time = 0.25;
                        data{5}.filter = 1;
                        data{5}.name = 'e';
                        
                        data{6}.pow = 700*ratio_power;
                        data{6}.acq = 5;
                        data{6}.time = 0.15;
                        data{6}.filter = 1;
                        data{6}.name = 'f';
                end
                
            case 'MBE2_4'
                switch detector
                    case 'vis'
                        data{1}.pow = 7.75*ratio_power ;
                        data{1}.acq = 50;
                        data{1}.time = 0.5;
                        data{1}.filter = NE30_exp_trans;
                        data{1}.name = 'a';
                        
                        data{2}.pow = 20.7*ratio_power ;
                        data{2}.acq = 50;
                        data{2}.time = 0.5;
                        data{2}.filter = NE40_exp_trans;
                        data{2}.name = 'b';
                        
                        data{3}.pow = 70.0*ratio_power ;
                        data{3}.acq = 50;
                        data{3}.time = 0.15;
                        data{3}.filter = NE40_exp_trans;
                        data{3}.name = 'c';
                        
                        data{4}.pow = 201*ratio_power;
                        data{4}.acq = 50;
                        data{4}.time = 0.3;
                        data{4}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{4}.name = 'd';
                        
                        data{5}.pow =450*ratio_power;
                        data{5}.acq = 50;
                        data{5}.time = 0.15;
                        data{5}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{5}.name = 'e';
                        
                    case 'IR'
                        data{1}.pow = 7.9*ratio_power ;
                        data{1}.acq = 5;
                        data{1}.time = 20;
                        data{1}.filter = 1;
                        data{1}.name = 'a';
                        
                        data{2}.pow = 20.8*ratio_power ;
                        data{2}.acq = 5;
                        data{2}.time = 5;
                        data{2}.filter = 1;
                        data{2}.name = 'b';
                        
                        data{3}.pow = 70.0*ratio_power ;
                        data{3}.acq = 5;
                        data{3}.time = 2;
                        data{3}.filter = 1;
                        data{3}.name = 'c';
                        
                        data{4}.pow = 200*ratio_power;
                        data{4}.acq = 5;
                        data{4}.time = 0.5;
                        data{4}.filter = 1;
                        data{4}.name = 'd';
                        
                        data{5}.pow =450*ratio_power;
                        data{5}.acq = 5;
                        data{5}.time = 0.2;
                        data{5}.filter = 1;
                        data{5}.name = 'e';
                end
                
            case 'MBE2_6'
                switch detector
                    case 'vis'
                        data{1}.pow = 7.2*ratio_power ;
                        data{1}.acq = 50;
                        data{1}.time = 0.3;
                        data{1}.filter = NE30_exp_trans;
                        data{1}.name = 'a';
                        
                        data{2}.pow = 19.8*ratio_power ;
                        data{2}.acq = 50;
                        data{2}.time = 0.3;
                        data{2}.filter = NE40_exp_trans;
                        data{2}.name = 'b';
                        
                        data{3}.pow = 69.1*ratio_power ;
                        data{3}.acq = 50;
                        data{3}.time = 0.5;
                        data{3}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{3}.name = 'c';
                        
                        data{4}.pow = 200*ratio_power;
                        data{4}.acq = 50;
                        data{4}.time = 0.5;
                        data{4}.filter = NE40_exp_trans.*NE20_exp_trans;
                        data{4}.name = 'd';
                        
                        data{5}.pow =450*ratio_power;
                        data{5}.acq = 50;
                        data{5}.time = 0.3;
                        data{5}.filter = NE40_exp_trans.*NE20_exp_trans;
                        data{5}.name = 'e';
                        
                        
                    case 'IR'
                        data{1}.pow = 7.9*ratio_power ;
                        data{1}.acq = 5;
                        data{1}.time = 20;
                        data{1}.filter = 1;
                        data{1}.name = 'a';
                        
                        data{2}.pow = 19.0*ratio_power ;
                        data{2}.acq = 5;
                        data{2}.time = 5;
                        data{2}.filter = 1;
                        data{2}.name = 'b';
                        
                        data{3}.pow = 69.0*ratio_power ;
                        data{3}.acq = 5;
                        data{3}.time = 1;
                        data{3}.filter = 1;
                        data{3}.name = 'c';
                        
                        data{4}.pow = 200*ratio_power;
                        data{4}.acq = 5;
                        data{4}.time = 0.3;
                        data{4}.filter = 1;
                        data{4}.name = 'd';
                        
                        data{5}.pow =450*ratio_power;
                        data{5}.acq = 5;
                        data{5}.time = 0.15;
                        data{5}.filter = 1;
                        data{5}.name = 'e';

                end
        end
        
    case '2019_12_10'
        switch sample
            case '1927_2_1'
                data{1}.pow = 9.21*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.5;
                data{1}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 25.8*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.15;
                data{2}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 50.0*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.3;
                data{3}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 101*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.15;
                data{4}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =200*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.5;
                data{5}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 301*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.5;
                data{6}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 500*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.5;
                data{7}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{7}.name = 'g';    
                
                data{8}.pow = 700*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.5;
                data{8}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{8}.name = 'h';    
                
                data{9}.pow = 900*ratio_power;
                data{9}.acq = 50;
                data{9}.time = 0.5;
                data{9}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{9}.name = 'i';
                
            case '1926_2_5'
                data{1}.pow = 9.37*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.3;
                data{1}.filter = NE40_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 26.1*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.5;
                data{2}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 49.7*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.3;
                data{3}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 100*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.15;
                data{4}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =200*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.3;
                data{5}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 301*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.15;
                data{6}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 400*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.5;
                data{7}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{7}.name = 'g';    

            case 'MBE1_2_1'
                data{1}.pow = 9.40*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.15;
                data{1}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 26.3*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.15;
                data{2}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 50.6*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.3;
                data{3}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 100*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.15;
                data{4}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =200*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.15;
                data{5}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 301*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.15;
                data{6}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 500*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.15;
                data{7}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{7}.name = 'g';    
                
                data{8}.pow = 700*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.15;
                data{8}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{8}.name = 'h';    
                
                data{9}.pow = 900*ratio_power;
                data{9}.acq = 50;
                data{9}.time = 0.15;
                data{9}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{9}.name = 'i';
        end
    case '2019_12_11'
        switch sample
            case '1873_1_4'
                data{1}.pow = 12*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.15;
                data{1}.filter = NE30_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 25*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.3;
                data{2}.filter = NE40_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 50*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.15;
                data{3}.filter = NE40_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 100*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.5;
                data{4}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =200*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.15;
                data{5}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 300*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.3;
                data{6}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 500*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.15;
                data{7}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{7}.name = 'g';    
                
                data{8}.pow = 700*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.3;
                data{8}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{8}.name = 'h';    
            
            case '1927_2_1'
            
                data{1}.pow = 13*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.15;
                data{1}.filter = NE30_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 25*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.3;
                data{2}.filter = NE40_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 50*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.15;
                data{3}.filter = NE40_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 100*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.5;
                data{4}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =200*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.15;
                data{5}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 300*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.3;
                data{6}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 500*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.15;
                data{7}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{7}.name = 'g';    
                
                data{8}.pow = 700*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.3;
                data{8}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{8}.name = 'h';    
                                
            case '1872_1_4'
            
                data{1}.pow = 12*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.5;
                data{1}.filter = NE10_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 25*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.15;
                data{2}.filter = NE10_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 50*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.3;
                data{3}.filter = NE20_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 100*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.3;
                data{4}.filter = NE30_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =200*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.15;
                data{5}.filter = NE30_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 300*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.5;
                data{6}.filter = NE40_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 500*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.3;
                data{7}.filter = NE40_exp_trans;
                data{7}.name = 'g';    
                
                data{8}.pow = 700*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.15;
                data{8}.filter = NE40_exp_trans;
                data{8}.name = 'h';    
                
            case '1926_2_5'
                data{1}.pow = 13*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.15;
                data{1}.filter = NE10_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 25*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.15;
                data{2}.filter = NE20_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 50*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.3;
                data{3}.filter = NE30_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 100*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.15;
                data{4}.filter = NE30_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =200*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.3;
                data{5}.filter = NE40_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 300*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.15;
                data{6}.filter = NE40_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 500*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.5;
                data{7}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{7}.name = 'g';    
                
                data{8}.pow = 700*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.5;
                data{8}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{8}.name = 'h';    
                
            case 'MBE1_2_1'
                
                data{1}.pow = 12*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.15;
                data{1}.filter = NE20_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 25*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.3;
                data{2}.filter = NE30_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 50*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.5;
                data{3}.filter = NE40_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 100*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.15;
                data{4}.filter = NE40_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =200*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.5;
                data{5}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 300*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.5;
                data{6}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 500*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.3;
                data{7}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{7}.name = 'g';    
                
                data{8}.pow = 700*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.15;
                data{8}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{8}.name = 'h';    
                
            case '1874_1_2'
            
                data{1}.pow = 12*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.3;
                data{1}.filter = NE20_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 25*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.15;
                data{2}.filter = NE20_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 50*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.15;
                data{3}.filter = NE30_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 100*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.5;
                data{4}.filter = NE40_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =200*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.15;
                data{5}.filter = NE40_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 300*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 1;
                data{6}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 500*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.5;
                data{7}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{7}.name = 'g';    
                
                data{8}.pow = 700*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.3;
                data{8}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{8}.name = 'h';                                   
        end
        
    case '2020_01_15'
        switch sample
            case 'MBE1_3_2'
                switch laser
                    case 532
                        
                        data{1}.pow = 0.78*ratio_power ;
                        data{1}.acq = 50;
                        data{1}.time = 0.3;
                        data{1}.filter = NE20_exp_trans;
                        data{1}.name = 'a';
                        
                        data{2}.pow = 2.45*ratio_power ;
                        data{2}.acq = 50;
                        data{2}.time = 0.15;
                        data{2}.filter = NE30_exp_trans;
                        data{2}.name = 'b';
                        
                        data{3}.pow = 7.74*ratio_power ;
                        data{3}.acq = 50;
                        data{3}.time = 0.15;
                        data{3}.filter = NE40_exp_trans;
                        data{3}.name = 'c';
                        
                        data{4}.pow = 25.5*ratio_power;
                        data{4}.acq = 50;
                        data{4}.time = 0.5;
                        data{4}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{4}.name = 'd';
                        
                        data{5}.pow =49.9*ratio_power;
                        data{5}.acq = 50;
                        data{5}.time = 0.3;
                        data{5}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{5}.name = 'e';
                        
                        data{6}.pow = 100*ratio_power;
                        data{6}.acq = 50;
                        data{6}.time = 0.15;
                        data{6}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{6}.name = 'f';
                        
                        data{7}.pow = 200*ratio_power;
                        data{7}.acq = 50;
                        data{7}.time = 0.15;
                        data{7}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{7}.name = 'g';
                        
                        data{8}.pow = 300*ratio_power;
                        data{8}.acq = 50;
                        data{8}.time = 0.3;
                        data{8}.filter = NE40_exp_trans.*NE20_exp_trans;
                        data{8}.name = 'h';
                        
                        data{9}.pow = 500*ratio_power;
                        data{9}.acq = 50;
                        data{9}.time = 0.15;
                        data{9}.filter = NE40_exp_trans.*NE20_exp_trans;
                        data{9}.name = 'i';
                        
                        data{10}.pow = 700*ratio_power;
                        data{10}.acq = 50;
                        data{10}.time = 0.15;
                        data{10}.filter = NE40_exp_trans.*NE20_exp_trans;
                        data{10}.name = 'j';
                        
                        data{11}.pow = 900*ratio_power;
                        data{11}.acq = 50;
                        data{11}.time = 0.15;
                        data{11}.filter = NE40_exp_trans.*NE20_exp_trans;
                        data{11}.name = 'k';
                        
                        data{12}.pow = 1100*ratio_power;
                        data{12}.acq = 50;
                        data{12}.time = 0.15;
                        data{12}.filter = NE40_exp_trans.*NE20_exp_trans;
                        data{12}.name = 'l';
                
                    case 638
                        
                        data{1}.pow = 12*ratio_power ;
                        data{1}.acq = 50;
                        data{1}.time = 0.15;
                        data{1}.filter = NE40_exp_trans;
                        data{1}.name = 'a';
                        
                        data{2}.pow = 25*ratio_power ;
                        data{2}.acq = 50;
                        data{2}.time = 0.5;
                        data{2}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{2}.name = 'b';
                        
                        data{3}.pow = 50*ratio_power ;
                        data{3}.acq = 50;
                        data{3}.time = 0.5;
                        data{3}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{3}.name = 'c';
                        
                        data{4}.pow = 100*ratio_power;
                        data{4}.acq = 50;
                        data{4}.time = 0.3;
                        data{4}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{4}.name = 'd';
                        
                        data{5}.pow =200*ratio_power;
                        data{5}.acq = 50;
                        data{5}.time = 0.15;
                        data{5}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{5}.name = 'e';
                        
                        data{6}.pow = 300*ratio_power;
                        data{6}.acq = 50;
                        data{6}.time = 0.15;
                        data{6}.filter = NE40_exp_trans.*NE10_exp_trans;
                        data{6}.name = 'f';
                        
                        data{7}.pow = 500*ratio_power;
                        data{7}.acq = 50;
                        data{7}.time = 0.3;
                        data{7}.filter = NE40_exp_trans.*NE20_exp_trans;
                        data{7}.name = 'g';
                        
                        data{8}.pow = 700*ratio_power;
                        data{8}.acq = 50;
                        data{8}.time = 0.3;
                        data{8}.filter = NE40_exp_trans.*NE20_exp_trans;
                        data{8}.name = 'h';
                end
                
            case 'MBE2_Q3_2_2'
                switch laser
                    case 532
                        
                        data{1}.pow = 2.51*ratio_power ;
                        data{1}.acq = 5;
                        data{1}.time = 30;
                        data{1}.filter = 1;
                        data{1}.name = 'a';
                        
                        data{2}.pow = 8.10*ratio_power ;
                        data{2}.acq = 10;
                        data{2}.time = 10;
                        data{2}.filter = 1;
                        data{2}.name = 'b';
                        
                        data{3}.pow = 22.5*ratio_power ;
                        data{3}.acq = 10;
                        data{3}.time = 5;
                        data{3}.filter = 1;
                        data{3}.name = 'c';
                        
                        data{4}.pow = 50.5*ratio_power;
                        data{4}.acq = 10;
                        data{4}.time = 2;
                        data{4}.filter = 1;
                        data{4}.name = 'd';
                        
                        data{5}.pow =100*ratio_power;
                        data{5}.acq = 10;
                        data{5}.time = 1;
                        data{5}.filter = 1;
                        data{5}.name = 'e';
                        
                        data{6}.pow = 200*ratio_power;
                        data{6}.acq = 10;
                        data{6}.time = 0.5;
                        data{6}.filter = 1;
                        data{6}.name = 'f';
                        
                        data{7}.pow = 300*ratio_power;
                        data{7}.acq = 10;
                        data{7}.time = 0.3;
                        data{7}.filter = 1;
                        data{7}.name = 'g';
                        
                        data{8}.pow = 500*ratio_power;
                        data{8}.acq = 10;
                        data{8}.time = 0.2;
                        data{8}.filter = 1;
                        data{8}.name = 'h';
                        
                        data{9}.pow = 700*ratio_power;
                        data{9}.acq = 10;
                        data{9}.time = 0.15;
                        data{9}.filter = 1;
                        data{9}.name = 'i';
                        
                        data{10}.pow = 900*ratio_power;
                        data{10}.acq = 10;
                        data{10}.time = 0.1;
                        data{10}.filter = 1;
                        data{10}.name = 'j';
                        
                        data{11}.pow = 1100*ratio_power;
                        data{11}.acq = 10;
                        data{11}.time = 0.1;
                        data{11}.filter = 1;
                        data{11}.name = 'k';
                        
                    case 638
                        
                        data{1}.pow = 12*ratio_power ;
                        data{1}.acq = 5;
                        data{1}.time = 30;
                        data{1}.filter = 1;
                        data{1}.name = 'a';
                        
                        data{2}.pow = 24*ratio_power ;
                        data{2}.acq = 5;
                        data{2}.time = 30;
                        data{2}.filter = 1;
                        data{2}.name = 'b';
                        
                        data{3}.pow = 50*ratio_power ;
                        data{3}.acq = 10;
                        data{3}.time = 15;
                        data{3}.filter = 1;
                        data{3}.name = 'c';
                        
                        data{4}.pow = 100*ratio_power;
                        data{4}.acq = 10;
                        data{4}.time = 5;
                        data{4}.filter = 1;
                        data{4}.name = 'd';
                        
                        data{5}.pow =200*ratio_power;
                        data{5}.acq = 10;
                        data{5}.time = 3;
                        data{5}.filter = 1;
                        data{5}.name = 'e';
                        
                        data{6}.pow = 300*ratio_power;
                        data{6}.acq = 10;
                        data{6}.time = 2;
                        data{6}.filter = 1;
                        data{6}.name = 'f';
                        
                        data{7}.pow = 500*ratio_power;
                        data{7}.acq = 10;
                        data{7}.time = 1;
                        data{7}.filter = 1;
                        data{7}.name = 'g';
                        
                        data{8}.pow = 700*ratio_power;
                        data{8}.acq = 10;
                        data{8}.time = 1;
                        data{8}.filter = 1;
                        data{8}.name = 'h';
                        
                end
        end
        
    case '2020_06_18'
        switch sample
            case '1925_4_2'
                data{1}.pow = 1.04*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.5;
                data{1}.filter = NE20_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 3.07*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.15;
                data{2}.filter = NE20_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 10.2*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.15;
                data{3}.filter = NE30_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 24.9*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.3;
                data{4}.filter = NE40_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =49.9*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.15;
                data{5}.filter = NE40_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 99.3*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.5;
                data{6}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 201*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.3;
                data{7}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{7}.name = 'g';
                
                data{8}.pow = 300*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.15;
                data{8}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{8}.name = 'h';
                
                data{9}.pow = 500*ratio_power;
                data{9}.acq = 50;
                data{9}.time = 0.15;
                data{9}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{9}.name = 'i';
                
                data{10}.pow = 700*ratio_power;
                data{10}.acq = 50;
                data{10}.time = 0.3;
                data{10}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{10}.name = 'j';
                
                data{11}.pow = 900*ratio_power;
                data{11}.acq = 50;
                data{11}.time = 0.15;
                data{11}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{11}.name = 'k';
                
                data{12}.pow = 1100*ratio_power;
                data{12}.acq = 50;
                data{12}.time = 0.15;
                data{12}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{12}.name = 'l';
                
            case 'MBE1_4_1'
                
                data{1}.pow = 1.07*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.5;
                data{1}.filter = NE30_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 3.20*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.5;
                data{2}.filter = NE40_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 10.3*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.15;
                data{3}.filter = NE40_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 25.2*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 1;
                data{4}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =49.9*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 1;
                data{5}.filter = NE40_exp_trans.*NE10_exp_trans;
                data{5}.name = 'e';
                
            case 'MBE3_Q2_1_1'
                
                data{1}.pow = 1.06*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.1;
                data{1}.filter = NE40_exp_trans;
                data{1}.name = 'a';
                
                data{2}.pow = 3.15*ratio_power ;
                data{2}.acq = 50;
                data{2}.time = 0.5;
                data{2}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{2}.name = 'b';
                
                data{3}.pow = 10.4*ratio_power ;
                data{3}.acq = 50;
                data{3}.time = 0.3;
                data{3}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{3}.name = 'c';
                
                data{4}.pow = 26.0*ratio_power;
                data{4}.acq = 50;
                data{4}.time = 0.1;
                data{4}.filter = NE40_exp_trans.*NE20_exp_trans;
                data{4}.name = 'd';
                
                data{5}.pow =51.3*ratio_power;
                data{5}.acq = 50;
                data{5}.time = 0.3;
                data{5}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{5}.name = 'e';
                
                data{6}.pow = 100.5*ratio_power;
                data{6}.acq = 50;
                data{6}.time = 0.2;
                data{6}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{6}.name = 'f';
                
                data{7}.pow = 200*ratio_power;
                data{7}.acq = 50;
                data{7}.time = 0.1;
                data{7}.filter = NE40_exp_trans.*NE30_exp_trans;
                data{7}.name = 'g';
                
                data{8}.pow = 300*ratio_power;
                data{8}.acq = 50;
                data{8}.time = 0.5;
                data{8}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{8}.name = 'h';
                
                data{9}.pow = 500*ratio_power;
                data{9}.acq = 50;
                data{9}.time = 0.5;
                data{9}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{9}.name = 'i';
                
                data{10}.pow = 700*ratio_power;
                data{10}.acq = 50;
                data{10}.time = 0.5;
                data{10}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{10}.name = 'j';
                
                data{11}.pow = 900*ratio_power;
                data{11}.acq = 50;
                data{11}.time = 0.5;
                data{11}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{11}.name = 'k';
                
                data{12}.pow = 1100*ratio_power;
                data{12}.acq = 50;
                data{12}.time = 0.5;
                data{12}.filter = NE40_exp_trans.*NE30_exp_trans.*NE10_exp_trans;
                data{12}.name = 'l';
                
            case 'MBE4_Q2_1_1'
                
                data{1}.pow = 1.06*ratio_power ;
                data{1}.acq = 50;
                data{1}.time = 0.3;
                data{1}.filter = NE30_exp_trans;
                data{1}.name = 'a';
        end
                
end

%% Data extraction
for index=1:length(data)
    for condition = ["bright", "dark"]
        %Finds the file
        if strcmp(detector,'IR')
            if laser==638
                filename = strcat(direction_name, sample, '_r_IR_', data{index}.name);
            else
                filename = strcat(direction_name, sample, '_IR_', data{index}.name);
            end
        else
            if laser==638
                 filename = strcat(direction_name, sample, '_r_', data{index}.name);
            else
                filename = strcat(direction_name, sample, '_', data{index}.name);
            end
        end
        if strcmp(condition, 'dark') == 1
            filename = strcat(filename,'_dark');
        end
        
        %Processes the file
        filename = strcat(filename,'.txt');
        fileID = fopen(filename,'r');
        if strcmp(detector,'IR')
            dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        else
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
        end
        fclose('all');
        lambda = dataArray{:,1}; % in nm
        
        %Which object to fill with the data
        if strcmp(condition, 'dark') == 1
            Int_dark{index} = zeros(length(lambda),1);
            for j=2:size(dataArray,2)-1
                Int_dark{index}(:) = Int_dark{index}(:)+dataArray{:,j};
            end
        else
            Int{index} = zeros(length(lambda),1);
            for j=2:size(dataArray,2)-1
                Int{index}(:) = Int{index}(:)+dataArray{:,j};
            end
        end
    end
end

%% Calibrated intensity

E=flipud(h*c./(q*lambda*1e-9)); %in eV

for index=1:length(data)
    power(index)=data{index}.pow;
    if strcmp(detector,'IR')
        Int_cal{index}=flipud(Int{index}-Int_dark{index})/(data{index}.acq*data{index}.time)./data{index}.filter;
    else
        Int_cal{index}=flipud(Int{index}-Int_dark{index}).*calib/(data{index}.acq*data{index}.time)./data{index}.filter;
    end
    Int_cal{index}=filloutliers(abs(Int_cal{index}),'linear','movmedian',5,'ThresholdFactor',0.5);
end

