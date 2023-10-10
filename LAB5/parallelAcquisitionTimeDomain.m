function CAF_and_corr = parallelAcquisitionTimeDomain(receivedSignal,codeReplica,samplingFreq,IntermediateFrequency,Tcoh)
    
    f_max = 5000; %Hz
    f_min = -f_max;
    f_delta = 2/(3*Tcoh);

    if length(codeReplica) < (1050*16)*Tcoh/1e-3 %GPS
        carrierLength = Tcoh; 
        L = Tcoh*samplingFreq;

    else  %GAL
        carrierLength = 4*Tcoh; 
        L = Tcoh*samplingFreq*4;

    end

    corr = zeros(round(2*f_max/f_delta)+1,L);

        for dop = 1: round(2*f_max/f_delta)+1
            IntermediateFrequencyR = IntermediateFrequency + (f_min + (dop-1)*f_delta);
            carrier = carrierOut(samplingFreq,IntermediateFrequencyR,carrierLength)';
            signal1 = receivedSignal.*carrier;
            signal1f = fft(signal1);
            signal2f = conj(fft(codeReplica));
            corr(dop,:) = ifft(signal1f.*signal2f)';
        end
        CAF = abs(corr).^2;
    CAF_and_corr = {CAF,corr};
end
