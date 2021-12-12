function [c, cost] = bp_salsa(y, A, AH, p, mu, Nit, lambda)

% c = BP(y, A, AH, p, mu, Nit)
%
% BASIS PURSUIT
% Minimize || c ||_1 such that y = A c
% where A AH = p I
%
% INPUT
%   A, AH - function handles for A and A^H
%   mu - Augmented Lagrangian parameter
%   Nit - Number of iterations
%
% OUTPUT
%   c : minimizing vector
%
% Use [c, cost] = BP(...) to obtain cost function per iteration.
% Use BP(..., lambda) to minimize || lambda .* c ||_1.

% Reference
% M. V. Afonso, J. M. Bioucas-Dias, and M. A. T. Figueiredo.
% Fast image recovery using variable splitting and constrained optimization.
% IEEE Trans. Image Process., 19(9):2345Ð2356, September 2010.

if nargin < 7
    lambda = 1;
end

% Initialization
c = AH(y);
d = zeros(size(c));
cost = zeros(1, Nit);

for i = 1:Nit
    u = soft(c + d, lambda/mu) - d;
    d = (1/p) * AH(y - A(u));
    c = d + u;
    cost(i) = sum(lambda(:) .* abs(c(:)));
end

