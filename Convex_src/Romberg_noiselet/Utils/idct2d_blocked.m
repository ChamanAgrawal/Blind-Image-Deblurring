% idct2d_blocked.m
%
% Inverse of dct2d_blocked.
%
% Written by: Justin Romberg, Georgia Tech
% Created: March 2007

function X = idct2d_blocked(D, blocksize)

[nr,nc] = size(D);

if ((mod(nr,blocksize) ~= 0) || (mod(nc,blocksize)~=0))
  error('blocksize must divide both height and width evenly');
end

br = nr/blocksize;
bc = nc/blocksize;

X = reshape(idct(reshape(reshape(idct(reshape(D, blocksize, br*nc)), nr, nc)', blocksize, bc*nr)), nr, nc)'; 

