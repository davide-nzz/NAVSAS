clc
close all
clear
addpath ('dataset/dataset')
codeReceived = load ('codeReceived.mat');
codeReceived = codeReceived.codeReceived(1,:)';

%% step1
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
signalLength = 1e-3; %ms
chipRate = 1.023e6;
carrier = carrierOut(samplingFrequency,intermediateFrequency,signalLength)';

%% step2
PRN = 3;
codeIn = GPSCode(PRN);
codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);
codeReplicaIF = codeReplica.*carrier; 

figure %time domain
subplot(4,1,1),plot(1:200,real(carrier(101:300)),'o-'),title('carrier for downconversion')
subplot(4,1,2),plot(1:200,codeReplica(101:300),'o-'),title('local code')
subplot(4,1,3),plot(1:200,real(codeReplicaIF(101:300)),'o-'),title('code at IF')
subplot(4,1,4),plot(1:200,codeReceived(101:300),'o-'),title('received signal')

figure %PSD
subplot(4,1,1)
[pxx,f] = pwelch(real(carrier),[],[],[],samplingFrequency);
plot(f,10*log10(pxx)),grid on,title('carrier for downconversion'),xlabel('Frequency (Hz)')
subplot(4,1,2)
[pxx,f] = pwelch(real(codeReplica),[],[],[],samplingFrequency);
plot(f,10*log10(pxx)),grid on,title('code'),xlabel('Frequency (Hz)')
subplot(4,1,3)
[pxx,f] = pwelch(real(codeReplicaIF),[],[],[],samplingFrequency);
plot(f,10*log10(pxx)),grid on,title('code at IF'),xlabel('Frequency (Hz)')
subplot(4,1,4)
[pxx,f] = pwelch(real(codeReceived),[],[],[],samplingFrequency);
plot(f,10*log10(pxx)),grid on,title('received signal'),xlabel('Frequency (Hz)')

%% step3
%circular correlation
for PRN = 3:3 % it's PRN 3
    codeIn = GPSCode(PRN);
    signal1 = generateLocalCode(codeIn, samplingFrequency, chipRate);
    signal2 = codeReceived;
    corr = circularCrossCorrelation(signal1,signal2);
    figure,plot(corr,'-'),title(['Circular correlation with PRN ' num2str(PRN) ' for cycle'])
end


%FFT based circular cross correlation
for PRN = 3:3
    codeIn = GPSCode(PRN);
    signal1 = generateLocalCode(codeIn, samplingFrequency, chipRate);
    signal2 = codeReceived;
    signal1f = fft(signal1);
    signal2f = conj(fft(signal2));
    corr = ifft(signal1f.*signal2f)/length(signal1);
    figure,plot(flip(corr),'-'),title(['Circular correlation with PRN ' num2str(PRN) ' fft-based'])
end


%% step4
PRN = 5;
codeIn = GPSCode(PRN);
samplingFrequency = 16.368e6; %Hz
intermediateFrequency = 4.092e6; %Hz
chipRate = 1.023e6;
Tcoh = 1e-3;
L = Tcoh*samplingFrequency;

fileName = 'SignalRX_1.bin';
[fid, message] = fopen(fileName, 'rb');
samplesToRead = L;
[receivedSignal, cntData] = fread(fid, samplesToRead, 'double');
fclose(fid);

codeReplica = generateLocalCode(codeIn, samplingFrequency, chipRate);

CAF = serialAcquisition(receivedSignal,codeReplica,samplingFrequency,intermediateFrequency,Tcoh);
figure
    f_max = 5000; %Hz, as inside serialAcquisition
    f_min = -f_max;
    f_delta = 2/(3*Tcoh);
surf(1:L,f_min:f_delta:f_max,CAF),title('Cross Ambiguty Function')
ylim([f_min,f_max]);
ylabel('Doppler Frequency Hz')
xlim([1 L]);
xlabel('Code Delay in Samples')

figure
plot(1:L,CAF(11,:))
title('1D CAF')
xlabel('Code Delay in Samples')

figure
plot(f_min:f_delta:f_max,CAF(:,9821))
title('1D CAF')
xlabel('Doppler Frequency Hz')

