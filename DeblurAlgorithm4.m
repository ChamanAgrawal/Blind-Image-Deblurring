function [xEst, wEst] = DeblurAlgorithm4(ww,xx,yy)
    global L B d0_h mu A y d0_h rho N
    L1 = size(xx,1); 
    L2 = size(xx,2); 
    L = L1*L2;
    
    [y,B,A,h,x,C,N,K,conv_wx_image,B_real] = find_subspace(ww,xx,yy);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    addpath('Non_convex_src/');
    [u,v] = nonConvexBD(ww,xx,yy,B,C,L,N,K,2000);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,l] = wavedec2(conv_wx_image,1,'db1');
    CC = @(x) waverec2(C*x,l,'db1');
    
    mat = @(x) reshape(x,L1,L2);
    BB = @(x) mat(B_real*x);
    u = abs(u);
    v = abs(v);
    xEst = CC(v);
    xEst = (x(1,1)/xEst(1,1))*(xEst-min(min(xEst)));  % Computing the estimate with an scaling factor
    wEst = BB(u);
end

function [y,B,A,h,x,C,N,K,conv_wx_image,B_real] = find_subspace(ww,xx,yy)
    L1 = size(xx,1); 
    L2 = size(xx,2); 
    L = L1*L2;
    
    % Helper functions
    mat = @(x) reshape(x,L1,L2);
    vec = @(x) x(:);

    % Computing Matrix B
    % Making matrix B with columns from identity matrix
    % Making vector h, so that w = Bh
    w_vec = vec(ww);
    Indw = zeros(L,1);
    Indw = abs(w_vec)>0;

    j = 1;
    K = sum(Indw);
    B_real = sparse(L,K);
    h = zeros(K,1);
    
    F = dftmtx(L)/sqrt(L);
    B = F(:,1:K);
    
    for i = 1:L
        if(Indw(i) == 1)
            B_real(i,j) = Indw(i);
            h(j) =w_vec(i);
            j = j+1;
        end
    end
    
    % Computing matrix C
    conv_wx_image = fftshift(mat(yy));
    [alpha_conv,~] = wavedec2(conv_wx_image,1,'db1');
    
    [alpha_x,~] = wavedec2(xx,1,'db1');
    
    alpha = alpha_x;
    Ind_alpha_conv = abs(alpha_conv)>0.0005*max(abs(alpha_conv));
    
    Ind_alpha_x = zeros(1,length(alpha)); % Do this if you want to kill support info. from original image
    Ind = ((Ind_alpha_conv>0)|(Ind_alpha_x>0)); % Taking union of both supports
    
    fprintf('Number of non-zeros in x estimated from the blurred image: %.3d\n', sum(Ind_alpha_conv));

    %%% Compute matrix C
    j = 1;
    N = sum(Ind);
    C = sparse(size(alpha,2),N); 
    for i = 1:size(alpha,2)
        if(Ind(i) == 1)
            C(i,j) = Ind(i);
            x(j) = alpha(i);
            j = j+1;
        end
    end
    x = x';

    h = h/norm(h);
    A = F*C;
    x = x/norm(x);
    y =diag(B*h*x'*A');
end
