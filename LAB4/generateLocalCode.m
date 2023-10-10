function y = generateLocalCode(codeIn, samplingFreq, chipRate)
    
    ratio = samplingFreq/chipRate;
    y = zeros(ratio*length(codeIn),1);
    N = length(codeIn);
    modulated = zeros(2*N,1);

    if length(codeIn) < 1050 % GPS

        for ii = 1:N
            y(1 + ratio*(ii-1) : ratio*ii) = codeIn(ii);
        end
    
    else % GAL

        for jj = 1:N %moltiplication with the subcarrier
            modulated(2*jj-1) = 1*codeIn(jj);
            modulated(2*jj) = -1*codeIn(jj);
        end

        for ii = 1:2*N
            y(1 + ratio/2*(ii-1) : ratio/2*ii) = modulated(ii);
        end  

    end
end