% bdct_linapprox_ordering.m
%
% Returns the index ordering for taking linear approximations using a
% blocked DCT.  DCT coefficients within a block are taken in the standard
% JPEG zigzag order (see jpgzzind.m).  "Extra" coefficients are assigned to
% DCT blocks block-columnwise (there is never a difference of more than 1
% terms taken between each of the blocks).
%
% Usage: order = bdct_linapprox_ordering(n, blocksize)
% n - sidelength (image is nxn)
% blocksize - sidelength of blocks
%
% Written by: Justin Romberg, Georgia Tech
% Created: April 2007
%

function order = bdct_linapprox_ordering(n, blocksize)

if (mod(n,blocksize) ~= 0)
  error('blocksize must divide height and width evenly');
end

b = n/blocksize;

[zblock, pairs] = jpgzzind(blocksize, blocksize);

bigind = reshape(1:n^2,n,n);
order = zeros(n^2,1);
for ii = 1:blocksize^2    
  order((ii-1)*b^2+1:ii*b^2) = reshape(bigind(pairs(ii,1):blocksize:n,pairs(ii,2):blocksize:n), b^2, 1);
end
