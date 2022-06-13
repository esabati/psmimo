function [V0,V1,V2] = algorithm1(V0_old,a,h,tau,gamma,rho,Pc,sigma2,Pmax,N,K)
V1_old = randn(N,N);
V2_old = randn(N,N);
lambda = 0.1;
l = 0;

while l < 10
    l = l + 1;
    [V0,V1,V2,t,u] = problem12(V1_old,V2_old,V0_old,a,h,tau,gamma,rho,Pc,sigma2,Pmax,N,K,lambda); % solving problem 12
    lambda = u/t;
    if isnan(lambda)
        disp('out');
        V0 = V0_old;
        V1 = V1_old;
        V2 = V2_old;
        break
    end
    V1_old = V1;
    V2_old = V2;
    V0_old = V0;
end
end

