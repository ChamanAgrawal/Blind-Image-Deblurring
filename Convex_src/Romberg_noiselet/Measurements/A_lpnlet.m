% A_lpnlet.m
% 
% Combine lowpass dct2 and noiselets.
% Usage: y = A_lpnlet(x, n, OM1, OM2)
%
% Written by: Justin Romberg, Georgia Tech, jrom@ece.gatech.edu
% Created: May 2007
%

function y = A_lpnlet(x, n, OM1, OM2)

y = [A_dct2(x, n, OM1); A_noiselet(x, OM2)];
