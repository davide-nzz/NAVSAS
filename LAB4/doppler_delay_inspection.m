clc 
clear 
close all
addpath ('dataset/dataset')

%-- there are 50 portions in total
%-- lab5 code is used since it's faster
%-- check delay and doppler variation over the portions
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
chipRate = 1.023e6;
Tcoh = 1e-3;
portions = 50;
index = zeros(50,2);

for portion = 1:50%portions
    L = Tcoh*samplingFrequency;
    f_max = 5000; %Hz, as inside parallelAcquisitionTimeDomain
    f_min = -f_max;
    f_delta = 2/(3*Tcoh);
    
    fileName = 'SignalRX_1.bin';
    [fid, message] = fopen(fileName, 'rb');
    samplesToRead = L*portion;
    [receivedSignal, cntData] = fread(fid, samplesToRead, 'double');
    fclose(fid);

        codeIn = GPSCode(5);
        codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);
        receivedSignalL = receivedSignal(1+L*(portion-1):L*portion);
        CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignalL,codeReplica,samplingFrequency,intermediateFrequency,Tcoh);

        maximum = max(max(CAF_and_corr{1}));
        [dop,tau] = find(CAF_and_corr{1}==maximum);
        index(portion,1) = dop;
        index(portion,2) = tau;
%     figure
%     [X,Y] = meshgrid(1:L,f_min:f_delta:f_max);
%     surf(X,Y,CAF_and_corr{1}),title(['GPS Cross Ambiguty Function fft PRN = 5, portion = ', num2str(portion)])
%     ylim([f_min,f_max]);
%     ylabel('Doppler Frequency Hz')
%     xlim([1 L]);
%     xlabel('Code Delay in Samples')
end

dops = -5000 * ones(16,1) + f_delta * (0: round(2*f_max/f_delta))'; %doppler values
% most of them have dops(12) as doppler and 9822 delay