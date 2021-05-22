% A_dct2.m
%
% Subset of DCT2 coefficients, vectorized.
% Usage : y = A_dct2(x, n, OMEGA)
%
% Written by: Justin Romberg, Georgia Tech
% Created: May 2007

function y = A_dct2(x, n, OMEGA)

y2 = dct2(reshape(x,n,n));
y = y2(OMEGA);
