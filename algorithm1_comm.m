function [V,EE] = algorithm1_comm(K,N,M,a,h,tau,sigma2,gamma,rho,xi,Pc,Pmax,epsilon)
V_old = zeros(N,N,K + 1);
for k = 1:1+K
    v_rand = ones(N,1).*exp(1j*2*pi*rand(N,1))/sqrt(N)*sqrt(Pmax/(1+K));
    V_old(:,:,k) = v_rand*v_rand';
end

lambda = 1;
l = 0;

t_old = 0;

while(1)
    l = l + 1;
    [V,t,u] = problem12_comm(K,N,M,a,h,tau,sigma2,gamma,rho,Pc,Pmax,lambda,V_old);
    V_old = V;
    
    lambda = u/t;
    if abs(t - t_old) < epsilon
        break
    end
    t_old = t;
end

sumV = zeros(N,N);
for k = 1:K + 1
    sumV = sumV + V(:,:,k);
end

SE = (trace(sumV)/rho + Pc)*t;
EE = SE/(trace(sumV)/rho + Pc + SE*xi);
end