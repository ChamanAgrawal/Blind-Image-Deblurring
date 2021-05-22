% dct2d_blocked.m
%
% Takes a blocked 2D discrete cosine transform of an image.
%
% Usage: D = dct2d_blocked(X, blocksize)
% X - nr x nc image (nr/blocksize and nc/blocksize must be integers);
% blocksize - break X into blocksize x blocksize segments, take DCT of each
% D - nr x nc images containing blocked DCTs
%
% Written by: Justin Romberg, Georgia Tech
% Created: March 2007

function D = dct2d_blocked(X, blocksize)

[nr,nc] = size(X);

if ((mod(nr,blocksize) ~= 0) || (mod(nc,blocksize)~=0))
  error('blocksize must divide both height and width evenly');
end

br = nr/blocksize;
bc = nc/blocksize;

D = reshape(dct(reshape(reshape(dct(reshape(X, blocksize, br*nc)), nr, nc)', blocksize, bc*nr)), nr, nc)'; 
