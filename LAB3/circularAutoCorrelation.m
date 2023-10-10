function y = circularAutoCorrelation(signal)
    N = length(signal);
    y = zeros(N,1);
    for s = 0:N-1
        shifted = zeros(N,1);
        shifted(s+1:N) = signal(1:N-s);
        shifted(1:s) = signal(N-s+1:N);
        y(s+1) = (signal'*shifted)/N;
    end
end