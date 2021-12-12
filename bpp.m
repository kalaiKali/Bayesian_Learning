close all
clear

I = sqrt(-1);

printme = @(txt) print('-deps', sprintf('figures/Example_BP_%s',txt));




M = 100;
N = 256;
y = rand(M,1);
c = (1/N)*AT(y, M, N);
y2 = A(c, M, N);
recon_err = y2 - y;   % recon_err : reconstruction error

fprintf('Maximum reconstruction error : %g\n', max(recon_err))
M = 100;
m = (0:M-1)';

f1 = 10.5;
x = exp(I*2*pi*f1/M*m);             % x : signal

figure(1)
clf
subplot(2,1,1)
plot(m, real(x), m, imag(x),'--')
xlabel('Time (samples)')
title('Signal');
legend('Real part','Imaginary part')
box off
ylim1 = [-1.4 1.8];
ylim(ylim1)
fprintf('signal')
X = fft(x);                         % X : DFT of x
figure
subplot(2,1,1)
stem(m,abs(X),'marker','none')
xlabel('Frequency (DFT index)')
title('(A) Fourier coefficients (DFT)');
box off
fprintf('DFT')
N = 256;
X = (1/N)*AT(x,M,N);

figure
clf
subplot(2,1,1)
stem(abs(X),'marker','none')
title('(B) Fourier coefficients (least square solution)');
xlabel('Frequency (index)')
box off
xlim([0 N])
fprintf('LeastSquares')
err = x - A(X,M,N);
err_max = max(abs(err));
fprintf('Least squares: Maximum reconstruction error = %g\n', err_max);
% Define functions (Matlab function handles)
H = @(x) A(x,M,N);
HT = @(x) AT(x,M,N);

% Define algorithm parameters

p = N;                  % p : Parseval constant
Nit = 100;              % Nit : number of iterations
mu = 5;                 % mu : ADMM parameter

% Run basis pursuit algorithm
[c, cost] = bp_salsa(x, H, HT, p, mu, Nit);
err = x - A(c,M,N);
err_max = max(abs(err));
fprintf('Basis pursuit: Maximum reconstruction error = %g\n', err_max);
figure
clf
subplot(3,3,[1 2 4 5])
plot(cost)
title('Cost function history');
xlabel('Iteration')
it1 = 4;
del = cost(it1) - min(cost);
ylim([min(cost)-0.1*del cost(it1)])
xlim([0 Nit])
box off
fprintf('CostFunction')
figure
clf
subplot(2,1,1)
stem(0:N-1, abs(c), 'marker','none')
title('(C) Fourier coefficients (basis pursuit solution)');
xlabel('Frequency (index)')
box off
xlim([0 N])
fprintf('BasisPursuit')
