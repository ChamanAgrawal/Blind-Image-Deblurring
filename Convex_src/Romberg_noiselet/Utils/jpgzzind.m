% jpgzzind.m
%
% Produces jpeg's zigzagging indices.
% Usage: z = jpgzzind(nr, nc)
% nr - number of rows in the block
% nc - number of columns in the block
% z - nc*nrx1 vector of index permutation
%     i.e.  X(z), where X is an nrxnc image will produce a column vector of
%           the entries of X in the desired zigzagging order.
%
% Written by: Justin Romberg, Georgia Tech
% Created: March 2007

function [z, pairs, map] = jpgzzind(nr, nc)

% states: 1 = UR1, 2 = UR2, 3 = DL, 4 = LL1, 5 = LL2, 6 = DR
state = 1;
ir = 1;
ic = 1;

map = zeros(nr,nc);
pairs = zeros(nr*nc, 2);
z = zeros(nr*nc,1);
current = 1;

while (current <= nr*nc)
  
  map(ir, ic) = current;
  pairs(current,:) = [ir ic];
  z(current) = nr*(ic-1) + ir;
  
  switch state
    
    case 1
      if (ic < nc)
        ic = ic + 1;
      else
        ir = ir + 1;
      end
      state = 2;
      current = current + 1;
      
    case 2
      ir = ir + 1;
      ic = ic - 1;
      state = 3;
      current = current + 1;
      
    case 3
      if ((ic > 1) && (ir < nr))
        ir = ir + 1;
        ic = ic - 1;
        current = current + 1;
      else
        state = 4;
      end
      
    case 4
      if (ir < nr)
        ir = ir + 1;
      else
        ic = ic + 1;
      end
      state = 5;
      current = current + 1;
      
    case 5
      ir = ir - 1;
      ic = ic + 1;
      state = 6;
      current = current + 1;
      
    case 6
      if ((ir > 1) && (ic < nc))
        ir = ir - 1;
        ic = ic + 1;
        current = current + 1;
      else
        state = 1;
      end
      
  end
  
  
end