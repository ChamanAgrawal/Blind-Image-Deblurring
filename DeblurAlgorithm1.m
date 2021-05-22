function [xEst, wEst] = DeblurAlgorithm1(w,x,y)
    % Summary of this function goes here
    % Detailed explanation goes here
    %%% Input Details
    % Input variables are as defined in the paper
    % w for blur kernel, x for original image and y for blurred image
    
    % Add dependencies path
    addpath(fullfile('Convex_src/minFunc'));
    addpath(fullfile('Convex_src/minFunc_2012'));
    addpath(fullfile('Convex_src/minFunc','compiled'));
    addpath(fullfile('Convex_src/minFunc','mex'));
    addpath('Convex_src/');
    
    L1 = size(x,1); 
    L2 = size(x,2); 
    L = L1*L2;
    
    % Helper functions
    mat = @(x) reshape(x,L1,L2);
    vec = @(x) x(:);

    % Computing Matrix B
    % Making matrix B with columns from identity matrix
    % Making vector h, so that w = Bh
    w_vec = vec(w);
    Indw = zeros(L,1);
    Indw = abs(w_vec)>0;
    figure;
    plot(Indw);
    title('Vertorised Blur Kernel non-zero coefficients');
    j = 1;
    K = sum(Indw);
    B = sparse(L,K);
    h = zeros(K,1);
    for i = 1:L
        if(Indw(i) == 1)
            B(i,j) = Indw(i);
            h(j) =w_vec(i);
            j = j+1;
        end
    end

    % Define function BB
    BB = @(x) mat(B*x);
    BBT = @(x) B'*vec(x);
    
    % Computing matrix C
    conv_wx_image = fftshift(mat(y));
    [alpha_conv,~] = wavedec2(conv_wx_image,4,'db1');
    figure;
    plot(alpha_conv); title('wavelet coefficients of the convolved image');
    
    [alpha_x,~] = wavedec2(x,4,'db1');
    
    alpha = alpha_x;
%     Ind = zeros(1,length(alpha));
    Ind_alpha_conv = abs(alpha_conv)>0.0005*max(abs(alpha_conv));
    % Ind_alpha_conv is support recovered from blurred image; For actual recovery without oracle
    % info
%     Ind_alpha_x = abs(alpha_x)>0.0005*max(abs(alpha_x));
    % Ind_alpha_x is support recovered from original image; For oracle assisted
    % recovery
    
    Ind_alpha_x = zeros(1,length(alpha)); % Do this if you want to kill support info. from original image
    % Ind_alpha_conv = zeros(1,length(alpha)); % Do this if you want to kill support info. from blurred image
    Ind = ((Ind_alpha_conv>0)|(Ind_alpha_x>0)); % Taking union of both supports
    
    fprintf('Number of non-zeros in x estimated from the blurred image: %.3d\n', sum(Ind_alpha_conv));

    %%% Compute matrix C
    j = 1;
    N = sum(Ind);
    C = sparse(size(alpha,2),N); 
    for i = 1:size(alpha,2)
        if(Ind(i) == 1)
            C(i,j) = Ind(i);
            m(j) = alpha(i);
            j = j+1;
        end
    end
    m = m';
    
    %%% Define function CC
    [~,l] = wavedec2(conv_wx_image,4,'db1');
    CC = @(x) waverec2(C*x,l,'db1');
    CCT = @(x) (C'*(wavedec2(x,4,'db1'))');
    
    x_hat = waverec2(C*m,l,'db1');
    fprintf('Origonal image vs Wavelet approximated image: %.3e\n', norm(x-x_hat,'fro')/norm(x,'fro'));
    figure;
    imagesc(x_hat), title('Approximation of original image from few coeffs'), colormap(gray), colorbar;
    
    %%% Convex Program for deconvolution
    [M,H] = blindDeconvolve_implicit_2D(vec(y),CC,BB,4,CCT,BBT);
    
    [UM,SM,VM] = svd(M,'econ');
    [UH,SH,VH] = svd(H,'econ');
    [U2,S2,V2] = svd(SM*VM'*VH*SH);
    
    %%% Estimates of m and h and recovery errors
    mEst=sqrt(S2(1,1))*UM*U2(:,1);
    hEst=sqrt(S2(1,1))*UH*V2(:,1);
    
    %%% Estimates of x and w
    xEst = CC(mEst);
    xEst = (x(1,1)/xEst(1,1))*(xEst-min(min(xEst)));  % Computing the estimate with an scaling factor
    wEst = BB(hEst);
    
    %%% Move this to error calc part
    fprintf('xEst error: %.3e\n', norm(x-xEst,'fro')/norm(x,'fro'));
    fprintf('h error: %.3e\n', norm(h-(hEst*h(1)/hEst(1))/norm(h)));
end





