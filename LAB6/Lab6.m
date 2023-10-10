clc
close all
clear 

%% step1
%GPS
PRN = 1;
codeIn = GPSCode(PRN);

% % FILTER
% wn = 1/2;
% order = 4;
% [Bf, Af] = butter(order, wn, 'low');
% codeIn = filter(Bf, Af, codeIn);

samplingFrequency = 16.368; %MHz
chipRate = 1.023; %MHz
promptCode = generateLocalCode(codeIn, samplingFrequency, chipRate);
L = length(promptCode);
deltas = [1;0.5;0.25];

for k = 1 : 3
    delta = deltas(k);
    deltan = samplingFrequency/chipRate*delta;
    lateCode = [promptCode(deltan/2+1:L);promptCode(1:deltan/2)];
    earlyCode = [promptCode(L-deltan/2+1:L);promptCode(1:L-deltan/2)];
    refleCode = 0.4*[promptCode(L-4+1:L);promptCode(1:L-4)];

    promptCorrelation = circularAutoCorrelation(promptCode);
    promptCorrelationf = circularCrossCorrelation(promptCode,refleCode);
    lateCorrelation = circularCrossCorrelation(lateCode,promptCode);
    lateCorrelationf = circularCrossCorrelation(lateCode,refleCode);
    earlyCorrelation = circularCrossCorrelation(earlyCode,promptCode);
    earlyCorrelationf = circularCrossCorrelation(earlyCode,refleCode);

%     promptCorrelation = promptCorrelation + promptCorrelationf; % to add multipath
%     lateCorrelation = lateCorrelation + lateCorrelationf; % to add multipath
%     earlyCorrelation = earlyCorrelation + earlyCorrelationf;% to add multipath

    promptCorrelation = [promptCorrelation(L-32:L);promptCorrelation(1:32)];
    lateCorrelation = [lateCorrelation(L-32:L);lateCorrelation(1:32)];
    earlyCorrelation = [earlyCorrelation(L-32:L);earlyCorrelation(1:32)];
    
    figure,title(['GPS correlators with Delta = ',num2str(delta)])
    x = -2:1/16:2;
    subplot(4,1,1),hold on,grid on,plot(x,promptCorrelation,'o-'),plot(x,lateCorrelation,'o-'),plot(x,earlyCorrelation,'o-'),title(['GPS correlators with Delta = ',num2str(delta)])
    xlabel('Chips')
    legend('Prompt','Late','Early')
    hold off
    subplot(4,1,2),hold on,grid on,title(['GPS coherent S-curve with Delta = ',num2str(delta)])
    plot(x,(earlyCorrelation-lateCorrelation),'o-')
    hold off
    subplot(4,1,3),hold on,grid on,title(['GPS normalized S-curve with Delta = ',num2str(delta)])
    plot(x,(0.5*(earlyCorrelation-lateCorrelation)./(earlyCorrelation+lateCorrelation)),'o-')
    hold off
    subplot(4,1,4),hold on,grid on,title(['GPS power normalized S-curve with Delta = ',num2str(delta)])
    plot(x,(0.5*(earlyCorrelation.^2-lateCorrelation.^2)./(earlyCorrelation.^2+lateCorrelation.^2)),'o-')
    hold off
end

GALILEO
galCodes = load('GalileoCodes.mat');
PRN = 1;
codeIn = galCodes.GalE1b(PRN,:);

% % FILTER
% wn = 1/2;
% order = 4;
% [Bf, Af] = butter(order, wn, 'low');
% codeIn = filter(Bf, Af, codeIn);

samplingFrequency = 16.368; %MHz
chipRate = 1.023; %MHz
promptCode = generateLocalCode(codeIn, samplingFrequency, chipRate);
L = length(promptCode);
deltas = [1;0.5;0.25];

for k = 1 : 3
    delta = deltas(k);
    deltan = samplingFrequency/chipRate*delta;
    lateCode = [promptCode(deltan/2+1:L);promptCode(1:deltan/2)];
    earlyCode = [promptCode(L-deltan/2+1:L);promptCode(1:L-deltan/2)];
    refleCode = 0.4*[promptCode(L-4+1:L);promptCode(1:L-4)];

    promptCorrelation = circularAutoCorrelation(promptCode);
    promptCorrelationf = circularCrossCorrelation(promptCode,refleCode);
    lateCorrelation = circularCrossCorrelation(lateCode,promptCode);
    lateCorrelationf = circularCrossCorrelation(lateCode,refleCode);
    earlyCorrelation = circularCrossCorrelation(earlyCode,promptCode);
    earlyCorrelationf = circularCrossCorrelation(earlyCode,refleCode);

%     promptCorrelation = promptCorrelation + promptCorrelationf; % to add multipath
%     lateCorrelation = lateCorrelation + lateCorrelationf; % to add multipath
%     earlyCorrelation = earlyCorrelation + earlyCorrelationf;% to add multipath
    
    promptCorrelation = [promptCorrelation(L-32:L);promptCorrelation(1:32)];
    lateCorrelation = [lateCorrelation(L-32:L);lateCorrelation(1:32)];
    earlyCorrelation = [earlyCorrelation(L-32:L);earlyCorrelation(1:32)];
    
    figure,title(['GALIELO correlators with Delta = ',num2str(delta)])
    x = -2:1/16:2;
    subplot(4,1,1),hold on,grid on,plot(x,promptCorrelation,'o-'),plot(x,lateCorrelation,'o-'),plot(x,earlyCorrelation,'o-'),title(['GALILEO correlators with Delta = ',num2str(delta)])
    xlabel('Chips')
    legend('Prompt','Late','Early')
    hold off
    subplot(4,1,2),hold on,grid on,title(['GALILEO coherent S-curve with Delta = ',num2str(delta)])
    plot(x,(earlyCorrelation-lateCorrelation),'o-')
    hold off
    subplot(4,1,3),hold on,grid on,title(['GALILEO normalized S-curve with Delta = ',num2str(delta)])
    plot(x,(0.5*(earlyCorrelation-lateCorrelation)./(earlyCorrelation+lateCorrelation)),'o-')
    hold off
    subplot(4,1,4),hold on,grid on,title(['GALILEO power normalized S-curve with Delta = ',num2str(delta)])
    plot(x,(0.5*(earlyCorrelation.^2-lateCorrelation.^2)./(earlyCorrelation.^2+lateCorrelation.^2)),'o-')
    hold off
end

