function y = carrierOut(samplingFreq,IntermediateFrequency,signalLength)
    A = sqrt(2);
    t = 1/samplingFreq : 1/samplingFreq : signalLength;
    y = A*exp(1i*2*pi*IntermediateFrequency*t);
end
