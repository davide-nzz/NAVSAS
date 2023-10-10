function y = linearCrossCorrelation(signal1,signal2)
    N1 = length(signal1); %equal to N2
    y = zeros(N1+1,1);
    for s = 0:N1
        shifted = zeros(N1,1);
        shifted(s+1:N1) = signal1(1:N1-s);
        y(s+1) = (signal2'*shifted)/(N1);
    end

end