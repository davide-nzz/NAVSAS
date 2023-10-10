function y = circularCrossCorrelation(signal1,signal2)
    N1 = length(signal1); %equal to N2
    y = zeros(N1,1);
    for s = 1:N1
        shifted = zeros(N1,1);
        shifted(s+1:N1) = signal1(1:N1-s);
        shifted(1:s) = signal1(N1-s+1:N1);
        y(s) = (signal2'*shifted)/N1;
    end
end