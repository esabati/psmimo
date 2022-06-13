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
xi = -26; % dynamic power consumption coefficient [dBm/bps]
eps = 0.001; % algorithm convergence accuracy

lambda = 1; % carrier wavelength
d = lambda/2; % spacing between adjacent antennas

sensing_directions = [-54, -18, 18, 54]; % sensing directions [degrees]
phi = [-30, 30]; % angle of departure from the BS to the users [degrees]

%% linear parameters %%
Pmax_lin = 10^((Pmax - 30)/10);
Pc_lin = 10^((Pc - 30)/10);

sigma2_lin = 10.^((sigma2 - 30)/10);
tau_lin = 10.^(tau/10);
gamma_lin = 10.^((gamma - 30)/10);
xi_lin = 10^((xi - 30)/10);

%% signals %%
s0 = randn(N,1); % dedicated radar signal
s = randn(K,1); % data signals
s = s/norm(s);
v = ones(N,K); % matrix containing the beamforming vectors for all users
x = v*s + s0; % transmited signal

V0 = eye(N); % dedicated radar signal autocorrelation matrix

%% channel %%
h = zeros(N,K); % matrix containing the channels from the BS to all users
% h(1,:) = ones(1,K); % LOS channel

%% steering vector %% 
Ntheta = 181; % number of sensing directions
theta = linspace(-90,90,Ntheta); % sensing directions

a = zeros(N,Ntheta); % steering vector
for i = 1:Ntheta
    for j = 0:N - 1
        a(j + 1,i) = exp(1i*2*pi*(d/lambda)*j*sin(theta(i)));
    end
end

a_sens = zeros(N,M); % sensing steering vector
for m = 1:M
    for i = 0:N - 1
        a_sens(i + 1,m) = exp(1i*2*pi*(d/lambda)*j*sin(sensing_directions(m)));
    end
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
    SNR(k) = (abs(conj(h(:,k))'*v(:,k))^2)/(aux + conj(h(:,k))'*V0*h(:,k) + sigma2_lin(k)); 
end

R = sum(log2(1 + SNR)); % system communication throughput

%% test %%
[V0,V1,V2] = algorithm1(V0,a_sens,h,tau_lin,gamma_lin,rho,Pc_lin,sigma2_lin,Pmax_lin,N,K,M);

%% Beampattern gain %%
p = zeros(1,Ntheta); % beampattern gains
for i = 1:Ntheta
    aux = V1 + V2 + V0;
    p(i) = conj(a(:,i))'*aux*a(:,i);
end

%% figures %%
figure();
plot(theta,abs(p));
hold on;
grid on;
plot([theta(1) theta(end)],[gamma_lin(1) gamma_lin(end)],':k');
plot(sensing_directions,zeros(1,M),'xr','MarkerSize',10,'LineWidth',2);
plot(phi,zeros(1,K),'og','MarkerSize',10,'LineWidth',2);
hold off;
xlim([-90 90]);

