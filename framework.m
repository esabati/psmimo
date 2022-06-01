clear;
close all;

%% parameters %%
N = 16; % number of antennas at the BS
K = 2; % number of users
M = 4; % number of sensing directions

Pmax = 30; % maximum transmit power at the BS [dBm]
Pc = 25; % circuit power of the system [dBm]

sigma2 = -80*ones(1,K); % noise power at the users [dBm]
tau = 5*ones(1,K); % required minimum SINR at the users [dB]
gamma = 20*ones(1,M); % required radar beampattern gain [dBm]
rho = 0.35; % amplifier efficiency at the BS

lambda = 1; % carrier wavelength
d = lambda/2; % spacing between adjacent antennas

%% signals %%
s0 = randn(N,1); % dedicated radar signal
s = randn(K,1); % data signals
v = ones(N,K); % matrix containing the beamforming vectors for all users
x = v*s + s0; % transmited signal

%% autocorrelation matrix %%
[sc,lags] = xcorr(s,s,N - 1,'biased');
r = sc(N:end);
V0 = toeplitz(r,conj(r)); % dedicated radar signal autocorrelation matrix

%% channel %%
h = zeros(N,K); % matrix containing the channels from the BS to all users
h(1,:) = ones(1,K); % LOS channel

%% steering vector %% 
Ntheta = 181; % number of sensing directions
theta = linspace(-90,90,Ntheta); % sensing directions

a = zeros(N,Ntheta); % steering vector
for i = 1:Ntheta
    for j = 0:N - 1
        a(j + 1,i) = exp(1i*2*pi*(d/lambda)*j*sin(theta(i)));
    end
end

%% beampattern gain %%
p = zeros(1,Ntheta); % beampattern gains
for i = 1:Ntheta
    aux = 0;
    for j = 1:K
        aux = aux + v(:,j)*conj(v(:,j))';
    end
    p(i) = conj(a(:,i))'*aux*a(:,i);
end

%% SNR %%
SNR = zeros(1,K); % vector containing each user's SNR
for k = 1:K
    aux = 0;
    for i = 1:K
        if i ~= k
            aux = aux + abs(conj(h(:,k))'*v(:,i))^2;
        end
    end
    SNR(k) = (abs(conj(h(:,k))'*v(:,k))^2)/(aux + conj(h(:,k))'*V0*h(:,k) + sigma2(k)); 
end


        


