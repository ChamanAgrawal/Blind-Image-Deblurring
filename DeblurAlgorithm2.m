function [xEst, wEst] = DeblurAlgorithm2(w,x,y,lambda,M,N,iterations)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('Prida_src/');



% read in image
% f = imread(fullfile(dataPath,fName));
% f = imresize(f,0.5);

% parameter settings
params.MK = M;
params.NK = N;
params.niters = iterations ;

[xEst, wEst] = blind_deconv(y,lambda,params);
end

