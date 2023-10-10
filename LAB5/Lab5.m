clc
close all
clear
addpath 'dataset/Real signals [December 13, 2019]'
galCode = load ('GalileoCodes.mat');

%% STEP 1

%--parallel acquisition GPS
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
chipRate = 1.023e6;
Tcoh = 1e-3;
L = Tcoh*samplingFrequency;
f_max = 5000; %Hz, as inside parallelAcquisitionTimeDomain
f_min = -f_max;
f_delta = 2/(3*Tcoh);

fileName = 'SignalRX_1.bin';
[fid, message] = fopen(fileName, 'rb');
samplesToRead = L;
[receivedSignal, cntData] = fread(fid, samplesToRead, 'double');
fclose(fid);


PRN = 5;
codeIn = GPSCode(PRN);
codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);
%tic
CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignal,codeReplica,samplingFrequency,intermediateFrequency,Tcoh);
%toc

figure
[X,Y] = meshgrid(1:L,f_min:f_delta:f_max);
surf(X,Y,CAF_and_corr{1}),title(['GPS Cross Ambiguty Function PRN ',num2str(PRN)])
ylim([f_min,f_max]);
ylabel('Doppler Frequency Hz')
xlim([1 L]);
xlabel('Code Delay in Samples')
shading interp

%close all

%--parallel acquisition GALILEO
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
chipRate = 1.023e6;
Tcoh = 1e-3;
L = Tcoh*samplingFrequency*4; %x4 GALILEO
f_max = 5000; %Hz, as inside parallelAcquisitionTimeDomain
f_min = -f_max;
f_delta = 2/(3*Tcoh);

fileName = 'SignalRX_2.bin';
[fid, message] = fopen(fileName, 'rb');
samplesToRead = L;
[receivedSignal, cntData] = fread(fid, samplesToRead, 'double');
fclose(fid);

PRN = 2;
codeIn = galCode.GalE1b(PRN,:);
codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);
CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignal,codeReplica,samplingFrequency,intermediateFrequency,Tcoh);

figure
[X,Y] = meshgrid(1:L,f_min:f_delta:f_max);
surf(X,Y,CAF_and_corr{1}),title(['GALILEO Cross Ambiguty Function FFT',num2str(PRN)]),
ylim([f_min,f_max]);
ylabel('Doppler Frequency Hz')
xlim([1 L]);
xlabel('Code Delay in Samples')
shading interp


%% STEP 2

%--parallel acquisition noisy GPS
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
chipRate = 1.023e6;
Tcoh = 1e-3;
L = Tcoh*samplingFrequency
f_max = 5000; %Hz, as inside parallelAcquisitionTimeDomain
f_min = -f_max;
f_delta = 2/(3*Tcoh)

fileName = 'SignalRX_3.bin';
[fid, message] = fopen(fileName, 'rb');
samplesToRead = L;
[receivedSignal, cntData] = fread(fid, samplesToRead, 'double');
fclose(fid);

for PRN = 6 %6,18,21
    codeIn = GPSCode(PRN);
    codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);
    CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignal,codeReplica,samplingFrequency,intermediateFrequency,Tcoh);
    
    figure
    [X,Y] = meshgrid(1:L,f_min:f_delta:f_max);
    surf(X,Y,CAF_and_corr{1}),title(['Noisy GPS CAF FFT PRN = ',num2str(PRN)])
    ylim([f_min,f_max]);
    ylabel('Doppler Frequency Hz')
    xlim([1 L]);
    xlabel('Code Delay in Samples')
    shading interp
    
end

%--parallel acquisition noisy GPS coherent time extension
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
chipRate = 1.023e6;
Tcoh = 1e-3;

for PRN = 6 %6,18,21
    Nc = 3;
    Tcoh_te = Tcoh*Nc;
    L = Tcoh_te*samplingFrequency;
    f_max = 5000; %Hz, as inside parallelAcquisitionTimeDomain
    f_min = -f_max;
    f_delta = 2/(3*Tcoh_te);
    
    fileName = 'SignalRX_3.bin';
    [fid, message] = fopen(fileName, 'rb');
    samplesToRead = L;
    [receivedSignal, cntData] = fread(fid, samplesToRead, 'double');
    fclose(fid);

    corr_sum = zeros(round(2*f_max/f_delta)+1,round(L/Nc));

    codeIn = GPSCode(PRN);
    codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);
    codeReplica = repmat(codeReplica,Nc,1);
    CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignal,codeReplica,samplingFrequency,intermediateFrequency,Tcoh_te);
    corr = CAF_and_corr{2};

    for ii = 1:Nc
        corr_sum = corr_sum + corr(:,1+L/Nc*(ii-1):L/Nc*ii);
    end
    CAF = abs(corr_sum/Nc).^2; 

    figure
    [X,Y] = meshgrid(1:L/Nc,f_min:f_delta:f_max);
    surf(X,Y,CAF),title(['Noisy GPS CAF FFT PRN = ',num2str(PRN) ,' with coherent time extension Nc = ', num2str(Nc)])
    ylim([f_min,f_max]);
    ylabel('Doppler Frequency Hz')
    xlim([1 L/Nc]);
    xlabel('Code Delay in Samples')
    shading interp
    
end

%--parallel acquisition noisy GPS non-coherent time extension
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
chipRate = 1.023e6;
Tcoh = 1e-3;


for PRN = 6 %6,18,21
    Nnc = 3;
    L = Tcoh*samplingFrequency*Nnc;
    f_max = 5000; %Hz, as inside parallelAcquisitionTimeDomain
    f_min = -f_max;
    f_delta = 2/(3*Tcoh);
    
    fileName = 'SignalRX_3.bin';
    [fid, message] = fopen(fileName, 'rb');
    samplesToRead = L;
    [receivedSignal, cntData] = fread(fid, samplesToRead, 'double');
    fclose(fid);

    CAF_sum = zeros(round(2*f_max/f_delta)+1,L/Nnc);

    codeIn = GPSCode(PRN);
    codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);

    for ii = 1:Nnc
        receivedSignalL = receivedSignal(1+L/Nnc*(ii-1):L/Nnc*ii);
        CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignalL,codeReplica,samplingFrequency,intermediateFrequency,Tcoh);
        CAF_sum = CAF_sum + CAF_and_corr{1};
    end
    CAF = CAF_sum/Nnc;

    figure
    [X,Y] = meshgrid(1:L/Nnc,f_min:f_delta:f_max);
    surf(X,Y,CAF),title(['Noisy GPS CAF FFT PRN = ',num2str(PRN) ,' with non-coherent time extension Nc = ', num2str(Nnc)])
    ylim([f_min,f_max]);
    ylabel('Doppler Frequency Hz')
    xlim([1 L/Nnc]);
    xlabel('Code Delay in Samples')
    shading interp
    
end

%--parallel acquisition noisy GPS coherent and non-coherent time extension
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
chipRate = 1.023e6;
Tcoh = 1e-3;

for PRN = 6 %6,18,21
    Nc = 3;
    Nnc = 3;
    Tcoh_te = Tcoh*Nc;
    L = Tcoh_te*samplingFrequency;
    f_max = 5000; %Hz, as inside parallelAcquisitionTimeDomain
    f_min = -f_max;
    f_delta = 2/(3*Tcoh_te);
    
    fileName = 'SignalRX_3.bin';
    [fid, message] = fopen(fileName, 'rb');
    samplesToRead = L*Nnc;
    [receivedSignal, cntData] = fread(fid, samplesToRead, 'double');
    fclose(fid);

    codeIn = GPSCode(PRN);
    codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);
    codeReplica = repmat(codeReplica,Nc,1);
    CAF_sum = zeros(round(2*f_max/f_delta)+1,round(L/Nc));

    for jj = 1:Nnc
        corr_sum = zeros(round(2*f_max/f_delta)+1,round(L/Nc));
            CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignal(1+L*(jj-1):jj*L),codeReplica,samplingFrequency,intermediateFrequency,Tcoh_te);
            corr = CAF_and_corr{2};
            for ii = 1:Nc
                corr_sum = corr_sum + corr(:,1+L/Nc*(ii-1):L/Nc*ii);
            end
        CAF = abs(corr_sum/Nc).^2; 
        CAF_sum = CAF_sum + CAF;
    end

    CAF = CAF_sum/Nnc;
    figure
    [X,Y] = meshgrid(1:L/Nc,f_min:f_delta:f_max);
    surf(X,Y,CAF),title(['Noisy GPS CAF FFT PRN = ',num2str(PRN) ,' with coherent and non-coherent time extension Nc = ', num2str(Nc)])
    ylim([f_min,f_max]);
    ylabel('Doppler Frequency Hz')
    xlim([1 L/Nc]);
    xlabel('Code Delay in Samples')
    shading interp
    
end

%% STEP 3
close all
clear
clc
galCode = load ('GalileoCodes.mat');

%--parallel acquisition noisy GPS coherent and non-coherent time extension
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
chipRate = 1.023e6;
Tcoh = 1e-3;

for PRN = 1:32
    Nc = 5;
    Nnc = 5;
    Tcoh_te = Tcoh*Nc;
    L = Tcoh_te*samplingFrequency;
    f_max = 5000; %Hz, as inside parallelAcquisitionTimeDomain
    f_min = -f_max;
    f_delta = 2/(3*Tcoh_te);
    
    fileName = 'd04_roomI';
    [fid, message] = fopen(fileName, 'rb');
    samplesToRead = L*Nnc;
    [receivedSignal, cntData] = fread(fid, samplesToRead, 'int8');
    fclose(fid);

    codeIn = GPSCode(PRN);
    codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);
    codeReplica = repmat(codeReplica,Nc,1);
    CAF_sum = zeros(round(2*f_max/f_delta)+1,round(L/Nc));

    for jj = 1:Nnc
        corr_sum = zeros(round(2*f_max/f_delta)+1,round(L/Nc));
            CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignal(1+L*(jj-1):jj*L),codeReplica,samplingFrequency,intermediateFrequency,Tcoh_te);
            corr = CAF_and_corr{2};
            for ii = 1:Nc
                corr_sum = corr_sum + corr(:,1+L/Nc*(ii-1):L/Nc*ii);
            end
        CAF = abs(corr_sum/Nc).^2; 
        CAF_sum = CAF_sum + CAF;
    end

    CAF = CAF_sum/Nnc;
    figure
    [X,Y] = meshgrid(1:L/Nc,f_min:f_delta:f_max);
    surf(X,Y,CAF),title(['Noisy GPS Cross Ambiguty Function fft PRN = ',num2str(PRN) ,' with coherent time extension Nc = ', num2str(Nc)])
    ylim([f_min,f_max]);
    ylabel('Doppler Frequency Hz')
    xlim([1 L/Nc]);
    xlabel('Code Delay in Samples')
    shading interp
end

%--parallel acquisition GALILEO

for PRN = 1:1
    Nc = 3;
    Nnc = 3;
    Tcoh_te = Tcoh*Nc;
    L = Tcoh_te*samplingFrequency*4; %x4 GALILEO
    f_max = 5000; %Hz, as inside parallelAcquisitionTimeDomain
    f_min = -f_max;
    f_delta = 2/(3*Tcoh_te);
    
    fileName = 'd01_mixto';
    [fid, message] = fopen(fileName, 'rb');
    samplesToRead = L*Nnc;
    [receivedSignal, cntData] = fread(fid, samplesToRead, 'int8');
    fclose(fid);

    codeIn = galCode.GalE1b(PRN,:);
    codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);
    codeReplica = repmat(codeReplica,Nc,1);
    CAF_sum = zeros(round(2*f_max/f_delta)+1,round(L/Nc));

    for jj = 1:Nnc
        corr_sum = zeros(round(2*f_max/f_delta)+1,round(L/Nc));
            CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignal(1+L*(jj-1):jj*L),codeReplica,samplingFrequency,intermediateFrequency,Tcoh_te);
            corr = CAF_and_corr{2};
            for ii = 1:Nc
                corr_sum = corr_sum + corr(:,1+L/Nc*(ii-1):L/Nc*ii);
            end
        CAF = abs(corr_sum/Nc).^2; 
        CAF_sum = CAF_sum + CAF;
    end

    CAF = CAF_sum/Nnc;
    figure
    [X,Y] = meshgrid(1:L/Nc,f_min:f_delta:f_max);
    surf(X,Y,CAF),title(['Noisy GPS Cross Ambiguty Function fft PRN = ',num2str(PRN) ,' with coherent time extension Nc = ', num2str(Nc)])
    ylim([f_min,f_max]);
    ylabel('Doppler Frequency Hz')
    xlim([1 L/Nc]);
    xlabel('Code Delay in Samples')
    shading interp
end