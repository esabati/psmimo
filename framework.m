clear;
close all;

%% parameters %%
N = 16; % number of antennas at the BS
K = 2; % number of users
M = 4; % number of sensing directions

Pmax = 10^(30/10)*0.001; % maximum transmit power at the BS
Pc = 10^(25.6/10)*0.001; % circuit power of the system

sigma2 = 10^(-80/10)*0.001*1e10; % noise power at the users
tau = 10^(5/10); % required minimum SINR at the users
gamma = 10^(20/10)*0.001; % required radar beampattern gain
rho = 0.35; % amplifier efficiency at the BS
xi = 10^(-26/10)*0.001; % dynamic power consumption coefficient [dBm/bps]
eps = 0.001; % algorithm convergence accuracy

lambda = 1; % carrier wavelength
d = lambda/2; % spacing between adjacent antennas

sensing_directions = [-54, -18, 18, 54]*pi/180; % sensing directions [degrees]
phi = [-30, 30]*pi/180; % angle of departure from the BS to the users [degrees]

%% channel %%
h = zeros(N,1,K); % matrix containing the channels from the BS to all users
alpha = 10^(-99/10)*1e10; % channel attenuation

for k = 1:K
    h(:,:,k) = sqrt(alpha)*exp(1i*2*pi*(d/lambda).*[0:N - 1]*sin(phi(k)));
end
%% steering vector %%
Ntheta = 181; % number of sensing directions
theta = linspace(-90,90,Ntheta)*pi/180; % sensing directions

a = zeros(N,1,Ntheta); % steering vector
for i = 1:Ntheta
    a(:,:,i) = (1/sqrt(N))*exp(1i*2*pi*(d/lambda).*[0:N - 1]*sin(theta(i)));
end

a_sens = zeros(N,1,M); % sensing steering vector
for m = 1:M
    a_sens(:,:,m) = (1/sqrt(N))*exp(1i*2*pi*(d/lambda).*[0:N - 1]*sin(sensing_directions(m)));
end

%% test %%
[V] = algorithm1(K,N,M,a_sens,h,tau,sigma2,gamma,rho,Pc,Pmax,eps);

%% beampattern gain %%
p = zeros(1,Ntheta); % beampattern gains
aux = V(:,:,1) + V(:,:,2) + V(:,:,3); 
for i = 1:Ntheta
    p(i) = a(:,i)'*aux*a(:,i);
end

%% figures %%
figure();
plot(theta*180/pi,abs(p));
hold on;
grid on;
plot([theta(1)*180/pi theta(end)*180/pi],[gamma(1) gamma(end)],':k');
plot(sensing_directions*180/pi,zeros(1,M),'xr','MarkerSize',10,'LineWidth',2);
plot(phi*180/pi,zeros(1,K),'og','MarkerSize',10,'LineWidth',2);
hold off;
xlim([-90 90]);
legend('ISAC signal (N = 16)','\Gamma thresholds','Target directions','User directions'); 
xlabel('\theta (degrees)');
ylabel('Beampattern gain');
