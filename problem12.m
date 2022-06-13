function [V0,V1,V2,t,u] = problem12(V1_old,V2_old,V0_old,a,h,tau,gamma,rho,Pc,sigma2,Pmax,N,K,M,lambda)
a1 = log(conj(h(:,1))'*V0_old*h(:,1) + conj(h(:,1))'*V2_old*h(:,1) + sigma2(1))/log(2);
a2 = log(conj(h(:,2))'*V0_old*h(:,2) + conj(h(:,2))'*V1_old*h(:,2) + sigma2(2))/log(2);

%% solution for each k %%
cvx_begin quiet
    variable V0(N,N) semidefinite
    variable V1(N,N) semidefinite
    variable V2(N,N) semidefinite
    variable t
    variable u
    maximize(t)
    R1 = log(conj(h(:,1))'*V0*h(:,1) + conj(h(:,1))'*V1*h(:,1) + conj(h(:,1))'*V2*h(:,1) + sigma2(1))/log(2);
    R2 = log(conj(h(:,2))'*V0*h(:,2) + conj(h(:,2))'*V1*h(:,2) + conj(h(:,2))'*V2*h(:,2) + sigma2(2))/log(2);
    subject to
        u >= (1/rho)*(trace(V0) + trace(V1) + trace(V2)) + Pc % (7a)
        conj(h(:,1))'*V1*h(:,1) >= tau(1)*(sigma2(1) + conj(h(:,1))'*V0*h(:,1) + conj(h(:,1))'*V2*h(:,1)) % (7c) k = 1
        conj(h(:,2))'*V2*h(:,2) >= tau(2)*(sigma2(2) + conj(h(:,2))'*V0*h(:,2) + conj(h(:,2))'*V1*h(:,2)) % (7c) k = 2
        trace(V1) + trace(V2) + trace(V0) <= Pmax % (7d)
        for m = 1:M
            real(conj(a(:,m))'*(V0 + V1 + V2)*a(:,m)) >= gamma(m) % (7e)
        end
        t >= 0 % (7f)
        u >= 0 % (7f)
cvx_end
disp('in')
end
