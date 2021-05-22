function [hEst,xEst] = nonConvexBD(w,xt,y,B,C,L,N,K,niters)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    global A B d y rho mu nue L h x
    F = dftmtx(L)/sqrt(L);
    A = F*C;
    B = B;
    
    %% Algorithm-1
    %Step-1
    A_star = operator_Astar(y);
%     
%     fprintf('size(A) = %i %i\n',size(A));
%     fprintf('size(B) = %i %i\n',size(B));
%     fprintf('size(A_star) = %i %i\n',size(A_star));
    %Step-2
    [U,S,V] = svd(A_star);
    d = S(1,1);
    h_hat_0 = U(:,1);
    x_hat_0 = V(:,1);
%     fprintf('length(h_hat_0) = %i\n',length(h_hat_0));
%     fprintf('length(x_hat_0) = %i\n',length(x_hat_0));
    
    
    %rho = d^2 + 2*norm(e)^2;
    rho = d^2/100;
    mu = 6*sqrt(L/((K+N)*log(L)^2));
    C_L = d*(N*log(L) + ((rho*L)/((d^2)*(mu^2))));
    nue = 1.0/C_L;
    
    %Step-3 (Solving optimization problem)    
    cvx_begin
        variable z(length(h_hat_0))
        minimize( norm(z - sqrt(d)*h_hat_0) )
        subject to 
            sqrt(L)*norm(B*z,Inf) <= 2*sqrt(d)*mu;
    cvx_end
    
    u_0 = z;
    v_0 = sqrt(d)*x_hat_0;
    
    %Step-4
    %fprintf('u_0 = %.3d',u_0);
    %fprintf('v_0 = %.3d',v_0);
    
    %% Algorithm-2
    
    u_t = u_0;
    v_t = v_0;    
    step_size = nue;
    
    for i = 1:niters
       [dF_tilda_h,dF_tilda_x] = grad(u_t,v_t);
       u_t = u_t - step_size*dF_tilda_h;
       v_t = v_t - step_size*dF_tilda_x;
%        if norm(h*x' - u_t*v_t','fro')/norm(u_t*v_t','fro') < 0.02
%            break;
%        end
    end
    
    %u_t = u_t/ norm(u_t);
    %v_t = v_t/ norm(u_t);
%     fprintf('norm u = %d',norm(u_t));
%     fprintf('norm v = %d',norm(v_t));
    hEst = (u_t);
    xEst = (v_t);
    
end

function val = operator_Astar(yy)
    global A B
    val = zeros(size(B,2),size(A,2));
    AT = A';
    BT = B';
    for i=1:length(yy)
        
        val = val + yy(i) .* BT(:,i) * AT(:,i)';
    end
end

function val = operator_A(Z)
    global A B
    val = diag(B * Z * A');
end

function res = g0_p(x)
    res = 2*sqrt(max(x-1,0)^2);
end

function val = helper(h)
    global B d mu L
    B_star = B';
    val = zeros(size(B,2),1);
    for i =1 : L
        val = val + g0_p(L*abs(B_star(:,i)'*h)^2/(8*d*mu^2))*B_star(:,i)*B_star(:,i)'*h;
    end
end

function [dF_tilda_h,dF_tilda_x] = grad(h,x)
    global d y rho mu
    dFh = operator_Astar(operator_A(h*x')-y)*x;
    dFx = operator_Astar(operator_A(h*x')-y)'*h;
    dGh = (rho/2.0*d)*((g0_p(norm(h)^2/2/d))*h + length(y)/(4*mu*mu)*helper(h));
    dGx = (rho/2.0*d) *(g0_p(norm(x)^2/2/d))*x;
    dF_tilda_h = dFh + dGh;
    dF_tilda_x = dFx + dGx;
end

