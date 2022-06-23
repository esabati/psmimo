function [V,t,u] = problem12(K,N,M,a,h,tau,sigma2,gamma,rho,Pc,Pmax,lambda,V_old)
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
            a_k = log(quad_form(h(:,:,k),sumV_old - V_old(:,:,k + 1)) + sigma2)/log(2); 
            b_k = 1/(quad_form(h(:,:,k),sumV_old - V_old(:,:,k + 1)) + sigma2)/log(2); 
            R_underline(1,k) = log(quad_form(h(:,:,k),sumV) + sigma2)/log(2) - (a_k + b_k*quad_form(h(:,:,k),(sumV - V(:,:,k + 1)) - (sumV_old - V_old(:,:,k + 1))));
        end
        lambda/2*t^2 + 1/(2*lambda)*u^2 <= sum(R_underline);
        
        u >= trace(sumV)/rho + Pc
        
        for k = 1:K
            quad_form(h(:,:,k),V(:,:,k + 1)) >= tau*(quad_form(h(:,:,k),sumV - V(:,:,k + 1)) + sigma2);
        end
        
        trace(sumV) <= Pmax % (7d)
        
        for m = 1:M
            quad_form(a(:,:,m),sumV) >= gamma % (7e)
        end
cvx_end
end
