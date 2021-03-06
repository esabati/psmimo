function [V,t,u] = problem12_comm(K,N,M,a,h,tau,sigma2,gamma,rho,Pc,Pmax,lambda,V_old)
cvx_begin sdp quiet
    variable t
    variable u
    variable V(N,N,1 + K) hermitian
    maximize t
    subject to
        sumV = zeros(N,N);
        sumV_old = zeros(N,N);
        
        for k = 1:K + 1
            sumV = sumV + V(:,:,k);
            sumV_old = sumV_old + V_old(:,:,k);
            V(:,:,k) >= 0 
        end
        
        t >= 0
        u >= 0
        
        clear R_underline
        for k = 1:K
            a_k = log(quad_form(h(:,:,k),sumV_old - V_old(:,:,k)) + sigma2)/log(2); 
            b_k = 1/(quad_form(h(:,:,k),sumV_old - V_old(:,:,k)) + sigma2)/log(2); 
            R_underline(1,k) = log(quad_form(h(:,:,k),sumV) + sigma2)/log(2) - (a_k + b_k*quad_form(h(:,:,k),(sumV - V(:,:,k)) - (sumV_old - V_old(:,:,k))));
        end
        lambda/2*t^2 + 1/(2*lambda)*u^2 <= sum(R_underline);
        
        u >= trace(sumV)/rho + Pc
        
        for k = 1:K
            quad_form(h(:,:,k),V(:,:,k)) >= tau*(quad_form(h(:,:,k),sumV - V(:,:,k)) + sigma2);
        end
        
        trace(sumV) <= Pmax % (7d)
cvx_end
end
