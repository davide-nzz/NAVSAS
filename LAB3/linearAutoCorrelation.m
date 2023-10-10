function y = linearAutoCorrelation(signal)
    N = length(signal);
    y = zeros(N+1,1);
    for s = 0:N
        shifted = zeros(N,1);
        shifted(s+1:N) = signal(1:N-s);
        y(s+1) = (signal'*shifted)/(N);
    end
end