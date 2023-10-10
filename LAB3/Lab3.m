clc
close all
clear 

%% STEP 1

%generate 32 Gold sequences

m = 10;
CODEs = zeros(1023,32);

for PRN = 1:32
    CODEs(:,PRN) = GPSCode(PRN);
end

% verify Gcodes properties
plusOnes = zeros(32,1);
minusOnes = zeros(32,1);



for PRN =  1:32
    plusOnes(PRN) = 0;
    minusOnes(PRN) = 0;
    for k = 1 : 2^(m)-1
        if CODEs(k,PRN) == 1
            plusOnes(PRN) = plusOnes(PRN)+1;
        else
            minusOnes(PRN) = minusOnes(PRN)+1;
        end
    end
end

%% STEP2

%linear auto corr
for PRN = 1:1
    signal = CODEs(:,PRN);
    corr = linearAutoCorrelation(signal);
    figure, plot([flip(corr(2:end));corr],'-'), grid on, title(['Linear auto correlation GPS PRN ' num2str(PRN)]), xlabel('Delay(chips)')
end

%circular auto corr
for PRN = 1:1
    signal = CODEs(:,PRN);
    corr = circularAutoCorrelation(signal);
    figure, plot([corr(round(length(corr)/2)+1:end);corr(1:round(length(corr)/2))],'-'), grid on, title(['Circular auto correlation GPS PRN ' num2str(PRN)]), xlabel('Delay(chips)')
end

[mm, ind_GPS] = max(corr);
corr(ind_GPS) = -1;
ratio_dB_GPS = 10*log10(1/max(corr))

%linear cross corr
for PRN = 1:1
    signal1 = CODEs(:,PRN);
    signal2 = CODEs(:,PRN+1);
    corr = linearCrossCorrelation(signal1,signal2);
    figure, plot([flip(corr(2:end));corr],'-'), grid on, title(['Linear cross-correlation between GPS PRN ' num2str(PRN) ' and GPS PRN ' num2str(PRN+1)]), xlabel('Delay(chips)')
end

%circular cross corr
for PRN = 1:1
    signal1 = CODEs(:,PRN);
    signal2 = CODEs(:,PRN+1);
    corr = circularCrossCorrelation(signal1,signal2);
    figure, plot([corr(round(length(corr)/2)+1:end);corr(1:round(length(corr)/2))],'-'), grid on, title(['Circular cross-correlation between GPS PRN ' num2str(PRN) ' and GPS PRN ' num2str(PRN+1)]), xlabel('Delay(chips)')
end

%% STEP 3

Gal_codes = load ('GalileoCodes.mat');

%linear auto corr
for codeNo = 1:1
    signal = Gal_codes.GalE1b(codeNo,:)';
    corr = linearAutoCorrelation(signal);
    figure, plot([flip(corr(2:end));corr],'-'), grid on, title(['Linear auto correlation GALILEO code No. ' num2str(codeNo)]), xlabel('Delay(chips)')
end

%circular auto corr
for codeNo = 1:1
    signal = Gal_codes.GalE1b(codeNo,:)';
    corr = circularAutoCorrelation(signal);
    figure, plot([corr(round(length(corr)/2)+1:end);corr(1:round(length(corr)/2))],'-'), grid on, title(['Circular auto correlation GALILEO code No. ' num2str(codeNo)]), xlabel('Delay(chips)')
end

[mm, ind_GPS] = max(corr);
corr(ind_GPS) = -1;
ratio_dB_GAL = 10*log10(1/max(corr))

%linear cross corr
for codeN0 = 1:1
    signal1 = Gal_codes.GalE1b(codeNo,:)';
    signal2 = Gal_codes.GalE1b(codeNo+1,:)';
    corr = linearCrossCorrelation(signal1,signal2);
    figure, plot([flip(corr(2:end));corr],'-'), grid on, title(['Linear cross-correlation between GALILEO code No. ' num2str(codeNo) 'and GALILEO code No.' num2str(codeNo+1)]), xlabel('Delay(chips)')
end

%circular cross corr
for codeN0 = 1:1
    signal1 = Gal_codes.GalE1b(codeNo,:)';
    signal2 = Gal_codes.GalE1b(codeNo+1,:)';
    corr = circularCrossCorrelation(signal1,signal2);
    figure, plot([corr(round(length(corr)/2)+1:end);corr(1:round(length(corr)/2))],'-'), grid on, title(['Circular cross-correlation of GALILEO code No. ' num2str(codeNo) 'and GALILEO code No.' num2str(codeNo+1)]), xlabel('Delay(chips)')
end

%% STEP4        

codeIn = CODEs(:,1);
samplingFreq = 16.368; %MHz
chipRate = 1.023; %MHz
codeOut_GPS = generateLocalCode(codeIn, samplingFreq, chipRate);

codeIn = Gal_codes.GalE1b(1,:)';
samplingFreq = 16.368; %MHz
chipRate = 1.023; %MHz
codeOut_GAL = generateLocalCode(codeIn, samplingFreq, chipRate);

%% STEP 5

figure
plot(codeOut_GPS(500:700))
title('GPS Code in time domain'), grid on, xlabel('seconds')

figure
plot(codeOut_GAL(500:700))
title('GALILEO Code in time domain'), grid on, xlabel('seconds')


figure
subplot(2,1,1),pwelch(codeOut_GPS)
title('GPS Code in frequency domain')
subplot(2,1,2),pwelch(codeOut_GAL)
title('GALILEO Code in frequency domain')

signal = codeOut_GPS;
corr_GPS = circularAutoCorrelation(signal);
corr_GPS_shifted = [corr_GPS(round(length(corr_GPS)/2)+1:end);corr_GPS(1:round(length(corr_GPS)/2))];
[val_GPS,ind_GPS] = max(corr_GPS_shifted);
signal = codeOut_GAL;
corr_GAL = circularAutoCorrelation(signal);
corr_GAL_shifted = [corr_GAL(round(length(corr_GAL)/2)+1:end);corr_GAL(1:round(length(corr_GAL)/2))];
[val_GAL,ind_GAL] = max(corr_GAL_shifted);

figure
subplot(2,1,1),plot(corr_GPS_shifted(ind_GPS-32:ind_GPS+32),'-'),title('CIrcular auto correlation local GPS PNR1'), grid on, xlabel('Normalized chip')

subplot(2,1,2),plot(corr_GAL_shifted(ind_GAL-32:ind_GAL+32),'-'),title('CIrcular auto correlation local GALILEO No1'), grid on, xlabel('Normalized chip')

