function [V,EE] = algorithm_radar(K,N,M,a,h,V_user,tau,sigma2,rho,Pc,xi,Pmax)
B = zeros(N,N);
for m = 1:M
    B = B + N*N*a(:,:,m)*a(:,:,m)';
end
B = real(B);

cvx_begin sdp quiet
    variable V(N,N) hermitian
    maximize trace(V*B)
    subject to
        sumV = V;
        for k = 1:K
            sumV = sumV + V_user(:,:,k);
        end
        
        SNR = zeros(1,K);
        for k = 1:K
            SNR(k) = quad_form(h(:,:,k),V_user(:,:,k));
            SNR(k) >= tau*(quad_form(h(:,:,k),sumV - V_user(:,:,k)) + sigma2)
        end
        V >= 0
        trace(V) <= Pmax
cvx_end

R = 0;
for k = 1:K
    R = R + log2(1 + SNR(k)/(quad_form(h(:,:,k),sumV - V_user(:,:,k)) + sigma2));
end
EE = R/((1/rho)*trace(sumV) + Pc + xi*R);

end

