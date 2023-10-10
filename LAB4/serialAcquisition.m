function CAF = serialAcquisition(receivedSignal,codeReplica,samplingFreq,IntermediateFrequency,Tcoh)

    L = Tcoh*samplingFreq;
    times = Tcoh*1000;

    codeReplicaR = [];
    for k = 1:times
        codeReplicaR = [codeReplicaR;codeReplica];
    end

    f_max = 5000; %Hz
    f_min = -f_max;
    f_delta = 2/(3*Tcoh);
    CAF = zeros(round(2*f_max/f_delta),L);
    
    for fd = 1: round(2*f_max/f_delta)+1
        IntermediateFrequencyR = IntermediateFrequency + (f_min + (fd-1)*f_delta);
        for tau = 1:L
            carrier = carrierOut(samplingFreq,IntermediateFrequencyR,Tcoh)';
                shifted = zeros(L,1);
                shifted(tau+1:L) =  codeReplicaR(1:L-tau);
                shifted(1:tau) = codeReplicaR(L-tau+1:L);
                codeReplicaIF = shifted.*(carrier);
            %codeReplicaIF_Re = shifted.*real(carrier);
            %codeReplicaIF_Im = shifted.*imag(carrier);
            CAF(fd,tau) = (abs(sum(codeReplicaIF.*receivedSignal)))^2;
        end
    end
end
