function [V_new,t,u] = problem12(V,V0,a,h,tau,rho,Pc,sigma2,Pmax,N,M,lambda)

cvx_begin
    variable V_new(N,N) semidefinite
    variable t
    variable u
    maximize(t)
    subject to
        (lambda/2)*t^2 + (1/(2*lambda))*u^2 <= 
        u >= (1/rho)*
    
cvx_end
end

