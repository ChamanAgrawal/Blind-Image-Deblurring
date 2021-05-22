clear all;
% camera_recovery_noiselet.m
%
% Create recovery curves for the 256x256 cameraman image
% Three approaches: 
%  linear approximation using dct2  (psnr_dct)
%  TV recovery from lowpass dct2 coefficients (psnr_tvlp)
%  TV recovery from 1000 lowpass + noiselet coefficients (psnr_cs)
%
% Written by: Justin Romberg, Georgia Tech, jrom@ece.gatech.edu
% Created: May 2007
%

addpath ./Measurements
addpath ./Optimization
addpath ./Utils


X = double(imread('cameraman.tif'));
x = X(:);
n = size(X,1);
N = n*n;

% for repeatable experiments
load RANDOM_STATES
rand('state', rand_state);
randn('state', randn_state);
%rand_state = rand('state');
%randn_state = randn('state');

% for linear approximation
lporder = bdct_linapprox_ordering(n, n);
% low pass DCT2, K1 = number of lowpass coefficients
K1 = 1000;
OM1 = lporder(1:K1);

% for random projection, avoid mean
q = randperm(N)';  
% K2 = number of auxiliary measurements
%  (either noiselet or more dct2)
K2 = 20000;  


% measurement ensemble
OM2 = q(1:K2);
Phi = @(z) A_lpnlet(z, n, OM1, OM2);
Phit = @(z) At_lpnlet(z, n, OM1, OM2);
% for linear and tvlp approximations
OMlin = lporder(1:K1+K2);
Phi2 = @(z) A_dct2(z, n, OMlin);
Phi2t = @(z) At_dct2(z, n, OMlin);

% take measurements
y = Phi(x);
y2 = Phi2(x);
  
% min l2 for cs
PPt = @(z) Phi(Phit(z));
x0 = Phit(cgsolve(PPt, y, 1e-8, 200));
  
% linear reconstruction
xlin = Phi2t(y2);
Xlin = reshape(xlin, n, n);

% parameters for optimization code
lbtol = 918;     % 1e-3*tv(X);
mu = 5;
cgtol = 1e-8;
cgmaxiter = 800;

% lowpass tv recovery
epsilon2 = 1e-3*norm(y2);
t0 = clock;
xlptv = tvqc_logbarrier(xlin, Phi2, Phi2t, y2, epsilon2, lbtol, mu, cgtol, cgmaxiter);
t1 = clock;
disp(sprintf('Elapsed time = %8.2f seconds', etime(t1,t0)));
Xlptv = reshape(xlptv, n, n);
  
% cs recovery
epsilon = 1e-3*norm(y);
t0 = clock;
xp = tvqc_logbarrier(x0, Phi, Phit, y, epsilon, lbtol, mu, cgtol, cgmaxiter);
t1 = clock;
disp(sprintf('Elapsed time = %8.2f seconds', etime(t1,t0)));
Xp = reshape(xp, n, n);

 
disp(sprintf('\n'));
disp(sprintf('K = %d + %d = %d', K1, K2, K1+K2));
disp(sprintf('DCT PSNR = %5.2f', psnr(X,Xlin)));
disp(sprintf('LPTV PSNR = %5.2f', psnr(X,Xlptv)));
disp(sprintf('CS PSNR = %5.2f', psnr(X,Xp)));
disp(sprintf('\n\n'));
