% At_lpnlet.m
% 
% Combine lowpass dct2 and noiselets (adjoint).
% Usage: w = At_lpnlet(y, n, OM1, OM2)
%
% Written by: Justin Romberg, Georgia Tech, jrom@ece.gatech.edu
% Created: May 2007
%

function w = At_lpnlet(y, n, OM1, OM2)

K1 = length(OM1);
K2 = K1 + length(OM2);
w = At_dct2(y(1:K1), n, OM1) + At_noiselet(y(K1+1:K2), OM2, n*n);

