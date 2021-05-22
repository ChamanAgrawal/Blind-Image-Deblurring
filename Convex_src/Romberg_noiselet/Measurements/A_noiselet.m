% A_noiselet.m
%
% Takes noiselet measurements on the specified set.
%
% Written by: Justin Romberg, Georgia Tech, jrom@ece.gatech.edu
% Created: May 2007
%

function y = A_noiselet(x, OMEGA)

n = length(x);
w = realnoiselet(x)/sqrt(n);
y = w(OMEGA);

