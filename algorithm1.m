function [V_new] = algorithm1(V0,a,h,tau,rho,Pc,sigma2,Pmax,N,M,epsilon)
V = randn(N,N);
V_new = zeros(N,N);
lambda = 1;
l = 0;

while V < epsilon
   l = l + 1;
   [V_new,t,u] = problem12(V,V0,a,h,tau,rho,Pc,sigma2,Pmax,N,M,lambda);
   lambda = u/t;
end
end

