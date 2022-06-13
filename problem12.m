function [V0,V1,V2,t,u] = problem12(N,M,a,h,tau,sigma2,gamma,rho,Pc,Pmax,lambda,V0_old,V1_old,V2_old)
a1 = log(conj(h(:,1))'*V0_old*h(:,1) + conj(h(:,1))'*V2_old*h(:,1) + sigma2(1))/log(2); % auxiliary variable for R1
a2 = log(conj(h(:,2))'*V0_old*h(:,2) + conj(h(:,2))'*V1_old*h(:,2) + sigma2(2))/log(2); % auxiliary variable for R2

disp('in');
cvx_begin quiet
    variable t
    variable u
    variable V0(N,N) semidefinite
    variable V1(N,N) semidefinite
    variable V2(N,N) semidefinite
    
    R1_aux = log(conj(h(:,1))'*V0*h(:,1) + conj(h(:,1))'*V1*h(:,1) + conj(h(:,1))'*V2*h(:,1) + sigma2(1))/log(2); % pre-calculation for R1
    R2_aux = log(conj(h(:,2))'*V0*h(:,2) + conj(h(:,2))'*V1*h(:,2) + conj(h(:,2))'*V2*h(:,2) + sigma2(2))/log(2); % pre-calculation for R2
    
    R1 = R1_aux - (a1 + ((log(exp(1))/log(2))/(2^a1))*(conj(h(:,1))'*(V0 - V0_old)*h(:,1) + conj(h(:,1))'*(V2 - V2_old)*h(:,1)));
    R2 = R2_aux - (a2 + ((log(exp(1))/log(2))/(2^a2))*(conj(h(:,2))'*(V0 - V0_old)*h(:,2) + conj(h(:,2))'*(V1 - V1_old)*h(:,2)));
    
    maximize(t)
    
    subject to
        (lambda/2)*t^2 + 1/(2*lambda)*u^2 <= R1 + R2 % (12)
        
        u >= (1/rho)*(trace(V0) + trace(V1) + trace(V2)) + Pc % (7a)
        
        conj(h(:,1))'*V1*h(:,1) >= tau(1)*(conj(h(:,1))'*V0*h(:,1) + conj(h(:,1))'*V2*h(:,1) + sigma2(1)) % (7c) k = 1
        conj(h(:,2))'*V2*h(:,2) >= tau(2)*(conj(h(:,2))'*V0*h(:,2) + conj(h(:,2))'*V1*h(:,2) + sigma2(2)) % (7c) k = 2
        
        trace(V0) + trace(V1) + trace(V2) <= Pmax % (7d)
        
        for m = 1:M
            (a(:,m)'*(V0 + V1 + V2)*a(:,m)) >= gamma(m) % (7e)
        end
        
        t >= 0 % (7f)
        u >= 0 % (7f)
cvx_end
end
