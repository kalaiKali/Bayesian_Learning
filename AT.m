function c = AT(y, M, N)
    c = fft([y; zeros(N-M, 1)]);
end