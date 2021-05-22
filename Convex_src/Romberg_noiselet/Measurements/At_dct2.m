% At_dct2.m
%
% Adjoint of A_dct2, vectorized.
% Usage : w = At_dct2(y, n, OMEGA)
%
% Written by: Justin Romberg, Georgia Tech
% Created: May 2007

function w = At_dct2(y, n, OMEGA)

y2 = zeros(n);
y2(OMEGA) = y;
w = reshape(idct2(y2), n^2, 1);
