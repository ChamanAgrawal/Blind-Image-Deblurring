function [snr] = psnr(original, noisy)
% Calculates PSNR between original image and noisy version

[height,width]=size(original);

error = abs(original - noisy);

enorm = sum(sum(error.^2));

snr = 10*log10(255*255*height*width/enorm);
