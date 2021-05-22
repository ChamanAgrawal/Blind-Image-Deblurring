function [Z,H] = blindDeconvolve_implicit_2D(conv_zh,C1,C2,maxrank,C1T,C2T,ZInit,HInit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function [Z,H] = blindDeconvolve(conv_zh,C1,C2,maxrank,ZInit,HInit)
%
%   Written by Ali Ahmed <alikhan@gatech.edu>
%   Last Updated:  November 20, 2012
%   Modified from BlindDeconvolve file written by Ben Recht.

%  Tries to find matrices Z and H such that
%       Z is n1 x maxrank
%       H is n2 x maxrank
%
%       conv_zh = sum_k circular_conv(Z(:,k),H(:,k))
%
%   Uses the method of multipliers to solve the algorithm
%
%   minimize ||L||^2 + ||R||^2
%   subject to conv_zh = sum_k circular_conv(Z(:,k),H(:,k))
%
%   based on the approach described in [1].  The inner iterations descend
%   the extended Lagrangian using minFunc.
%
%   Inputs ZInit and HInit are optional.  They give starting points for the
%   algorithm.
%
%
%   References:
%
%   [1] Samuel Burer and Renato D. C. Monteiro.  "A nonlinear programming
%   algorithm for solving semidefinite programs via low-rank
%   factorization." Mathematical Programming (Series B), Vol. 95, 
%   2003. pp 329-357.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parameters are
%
%   maxOutIter = maximum number of outer iterations
%
%   rmseTol = the root mean square of the errors must drop below
%   rmseTol before termination
%
%   sigmaInit = starting value for the Augmented Lagrangian parameter
%
%   LR1 = the feasibility must drop below LR1 times the previous
%   feasibility before updating the Lagrange multiplier.
%   Otherwise we update sigma.
%
%   LR2 = the amount we update sigma by if feasibility is not improved
%
%   progTol = like LR1, but much smaller.  If feasibility is not
%   improved by progtol for numbaditers iterations, then we decide
%   that further progress cannot be made and quit
%
%   numbaditers = the number of iterations where trivial progress
%   is made improving the feasibility of L and R before algorithm
%   termination.    

        
pars.maxOutIter = 25;
pars.rmseTol = 1e-8;
pars.sigmaInit = 1e4;
pars.LR1 = 0.25;
pars.LR2 = 10;
pars.progTol = 1e-3;
pars.numbaditers = 6;

% options for minFunc
options = [];
options.display = 'none';
options.maxFunEvals = 50000;
options.Method = 'lbfgs';


%% Derived constants:
siglen = length(conv_zh);
p = siglen;
p1 = sqrt(p);
p2 = sqrt(p);
mat = @(x) reshape(x,p1,p2);
vec = @(x) x(:);
n1 = length(C1T(mat(conv_zh))); % There must be a better way to do this...
n2 = length(C2T(mat(conv_zh)));


% we actually run the equality constraint in the frequency domain:
meas_fft = vec(fft2(mat(conv_zh)))/siglen;

% If starting points are specified, use them. otherwise use the default.
% A better default might speed things up.  Not sure.
if nargin>6 & ~isempty(LInit),
Z = ZInit;
else
Z = 1e-2*randn(n1,maxrank);
end

if nargin>7 & ~isempty(RInit),
H = HInit;
else
H = 1e-2*randn(n2,maxrank);
end
  
% y is the Lagrange multiplier
y = zeros(siglen,1);
% sigma is the penalty parameter
sigma = pars.sigmaInit;

% compute initial infeasibility
dev = zeros(p,1);
for i = 1:maxrank
    dev = dev+vec(fft2(mat(C1(Z(:,i)))).*fft2(mat(C2(H(:,i)))));
end
dev = vec(ifft2(mat(dev))) - conv_zh;
vOld = norm(dev,'fro')^2;
  
v = vOld;
badcnt = 0;
T0 = clock;
  
fprintf('|      |          |          | iter  | tot   |\n');
fprintf('| iter |  rmse    |  sigma   | time  | time  |\n');
for j=1:46, fprintf('-'); end
fprintf('\n');
  
iterCount = 0;

for outIter=1:pars.maxOutIter,
    
    T1 = clock;

    % minimize the Augmented Lagrangian using BFGS
    [x,mval] = minFunc(@subproblem_cost,[Z(:);H(:)],options,C1,C2,C1T,C2T,maxrank,meas_fft,y,sigma,siglen,n1,n2);
    Z = reshape(x(1:(n1*maxrank)),n1,maxrank);
    H = reshape(x((n1*maxrank) + (1:(n2*maxrank))),n2,maxrank);
    
    % compute the equation error
    dev = zeros(p,1);
    for i = 1:maxrank
    dev = dev+vec(fft2(C1(Z(:,i))).*fft2(C2(H(:,i))))/siglen;
    end
    dev = dev - meas_fft;
    vLast = v;
    v = norm(dev)^2; % v is sum of the squares of the errors
    
    % if unable to improve feasibility for several iterations, quit.
    if abs(vLast-v)/vLast<pars.progTol,
      badcnt = badcnt+1;
      if badcnt>pars.numbaditers,
        fprintf('\nunable to make progress. terminating\n');
        break;
      end
    else
      badcnt = 0;
    end
    
    % print diagnostics
    fprintf('| %2d   | %.2e | %.2e |  %3.0f  |  %3.0f  |\n',...
            outIter,sqrt(v/siglen),sigma,etime(clock,T1),etime(clock,T0));

    % if solution is feasible to specified tolerance, we're done
    if sqrt(v/siglen)<pars.rmseTol,
      break;
    % if normed feasibility is greatly reduced, update Lagrange multipliers
    elseif v < pars.LR1*vOld
      y = y - sigma*dev;
      vOld = v;
    % if normed feasibility is not reduced, increase penalty term
    else
      sigma = pars.LR2*sigma;    
    end
    
  end
  % print final diagnostics
  fprintf('elapsed time: %.0f seconds\n',etime(clock,T0));
  
function [mval,g]=subproblem_cost(x,C1,C2,C1T,C2T,maxrank,meas_fft,y,sigma,siglen,n1,n2)
% This helper function computes the value and gradient of the augmented
% Lagrangian  
    
    
    
    p1 = sqrt(siglen);
    p2 = p1;
    
    Z = reshape(x(1:(n1*maxrank)),n1,maxrank);
    H = reshape(x((n1*maxrank) + (1:(n2*maxrank))),n2,maxrank);
    mat = @(x) reshape(x,p1,p2);
    vec = @(x) x(:);
    % compute equation error
    dev = zeros(siglen,1);
    for i = 1:maxrank
    dev = dev+vec(fft2(C1(Z(:,i))).*fft2(C2(H(:,i))))/siglen;
    end    
    dev = dev - meas_fft;
    % compute the cost of the extended Lagrangian penalty function
    mval = norm(Z,'fro')^2+norm(H,'fro')^2 - 2*real(y'*dev) + sigma*norm(dev,'fro')^2; 
    
    % compute the gradient
    yhat = y - sigma*dev;
    temp1 = zeros(siglen,maxrank);
    temp2 = zeros(siglen,maxrank);
    for i = 1:maxrank 
       temp1(:,i) = vec(fft2(C2(H(:,i))));
       temp2(:,i) = vec(ifft2(C1(Z(:,i))));
    end
    temp3 = conj((yhat)*ones(1,maxrank)).*temp1;
    temp4 = (yhat)*ones(1,maxrank).*temp2;
    temp5 = zeros(n1,maxrank);
    temp6 = zeros(n2,maxrank);
    for i = 1:maxrank
       temp5(:,i) =  C1T(fft2(mat(temp3(:,i))));
       temp6(:,i) = C2T(ifft2(mat(temp4(:,i))));
    end
    adjoint_times_H = temp5/siglen;
    adjoint_times_Z = temp6*siglen;
    
    GradZ = 2*(Z - adjoint_times_H); 
    GradH = 2*(H - adjoint_times_Z);
    
    g = real([GradZ(:);GradH(:)]);