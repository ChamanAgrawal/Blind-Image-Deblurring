function howtousenoiseletcode_v1
% uses original Romberg C code

addpath '..\noiselet'
N = 256;
disp(' ');
fprintf(1,'\nrunning C code version of noiselets, N = %d',N);

s = RandStream('mt19937ar','Seed', 1234);
RandStream.setDefaultStream(s);

X = rand (N);    % input image
X = X(:);
M = N^2/4;  % number of compressive measurements 
% for full recovery M = length of X (||X||), for CS, M << ||X||

%% Make Compressive Observations 
% CS vars. 
% impt: These must be within scope for all code in this function
[N,junk]= size(X);     % 65536 x #rangebins for 256 x256 scenery
q       = randperm(N)';    % makes column vector of random integers 1:N
OM2     = q(1:N);          % vector of random subset of integers 1:N
Phi     = @(z) A_noiselet (z, OM2);
Phit    = @(z) At_noiselet(z, OM2, N);  % N conveys size of recovered x vector
A       = @(z) Phi( z );   
At      = @(z) Phit(z );

Y      = A(X(:));  % vector operand for noiselets  .  Fails for matrix input.
Xhat   = At(Y(:) );

error = abs(X - Xhat);
enorm = norm(error)/norm( X(:) ); 
fprintf(1,'\nenorm = %g',enorm);
d = 1;
end

