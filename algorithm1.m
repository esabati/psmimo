function [V0,V1,V2] = algorithm1(N,M,a,h,tau,sigma2,gamma,rho,Pc,Pmax)
V0_old = zeros(N,N);
V1_old = zeros(N,N);
V2_old = zeros(N,N);
lambda = 1;
l = 0;

while l < 10
    l = l + 1;
    [V0,V1,V2,t,u] = problem12(N,M,a,h,tau,sigma2,gamma,rho,Pc,Pmax,lambda,V0_old,V1_old,V2_old);
    if isnan(u) || isnan(t)
        V0 = V0_old;
        V1 = V1_old;
        V2 = V2_old;
        disp('out');
        break
    end
    V0_old = V0;
    V1_old = V1;
    V2_old = V2;
    
    lambda = u/t;
end
end