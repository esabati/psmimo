clear all;
close all;

N = 16; % number of antennas at the BS
K = 2; % number of users
M = 4; % number of sensing directions
lambda = 1; % carrier wavelength
d = lambda/2; % spacing between adjacent antennas

s0 = zeros(N,1); % dedicated radar signal
s = round(randn(K,1)); % data signals
V = ones(N,K); % matrix containing the beamforming vectors for all users
x = V*s + s0; % transmited signal

V0 = toeplitz(xcorr(s0)); % dedicated radar signal autocorrelation matrix

h = zeros(N,K); % matrix containing the channels from the BS to all users

Ntheta = 100; % number of sensing directions
theta = linspace(-90,90,Ntheta); % sensing directions

a = zeros(N,Ntheta); % steering vector

for i = 1:Ntheta
    for j = 0:N - 1
        a(j + 1,i) = exp(1i*2*pi*(d/lambda)*j*sin(theta(i)));
    end
end

p = zeros(1,Ntheta); % beampattern gains

for i = 1:Ntheta
    aux = 0;
    for j = 1:K
        aux = aux + V(:,j)*conj(V(:,j))';
    end
    p(i) = conj(a(:,i))'*aux*a(:,i);
end

        


