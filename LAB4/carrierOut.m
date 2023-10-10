function y = carrierOut(samplingFreq,IntermediateFrequency,signalLength)
    A = sqrt(2);
    t = 0 : 1/samplingFreq : signalLength-1/samplingFreq;
    y = A*exp(1i*2*pi*IntermediateFrequency*t);
end
