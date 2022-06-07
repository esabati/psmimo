function [V_new,t,u] = problem12(V,V0,a,h,tau,rho,Pc,sigma2,Pmax,N,M,lambda)

cvx_begin sdp
    variable V_new semidefinite
    variable t
    variable u
cvx_end
end

