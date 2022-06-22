function [V] = algorithm_radar(N,M,a,Pmax)
B = zeros(N,N);
for m = 1:M
    B = B + N*a(:,:,m)*a(:,:,m)';
end
B = real(B);

cvx_begin sdp quiet
    variable V(N,N) hermitian
    maximize trace(V*B)
    subject to
        trace(V) == Pmax
        V >= 0
cvx_end
end

