clear all;
clc;
close all;

%% parameters %%
N = 16; % number of antennas at the BS
K = 2; % number of users
M = 2; % number of sensing directions

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

sensing_directions = [-30, 30]*pi/180; % sensing directions [degrees]
phi = [-30, 30]*pi/180; % angle of departure from the BS to the users [degrees]

%% channel %%
h = zeros(N,1,K); % matrix containing the channels from the BS to all users
alpha = 10^(-99/10)*1e10; % channel attenuation

for k = 1:K
    h(:,:,k) = sqrt(alpha)*exp(1i*2*pi*(d/lambda).*[0:N - 1]*sin(phi(k)));
end
%% steering vector %%
Ntheta = 1810; % number of sensing directions
theta = linspace(-90,90,Ntheta)*pi/180; % sensing directions

a = zeros(N,1,Ntheta); % steering vector
for i = 1:Ntheta
    a(:,:,i) = (1/sqrt(N))*exp(1i*2*pi*(d/lambda).*[0:N - 1]*sin(theta(i)));
end

a_sens = zeros(N,1,M); % sensing steering vector
for m = 1:M
    a_sens(:,:,m) = (1/sqrt(N))*exp(1i*2*pi*(d/lambda).*[0:N - 1]*sin(sensing_directions(m)));
end

%% algorithm call %%
gamma = 0.001*10.^([10 12 14 16 18 20]/10);
EE_opt = zeros(1,length(gamma));

for i = 1:length(gamma)
    disp(i)
    [V,EE] = algorithm1(K,N,M,a_sens,h,tau,sigma2,gamma(i),rho,xi,Pc,Pmax,eps);
    EE_opt(i) = EE;
end
[V_comm,EE_comm] = algorithm1_comm(K,N,M,a_sens,h,tau,sigma2,gamma(end),rho,xi,Pc,Pmax,eps);
[V_radar,EE_radar] = algorithm_radar(K,N,M,a_sens,h,V(:,:,2:K+1),tau,sigma2,rho,Pc,xi,Pmax);

for k = 2:K + 1
    V_radar(:,:,k) = V(:,:,k);
end
%% SNR %%

sumV = sum(V,3);
sumV_comm = sum(V_comm,3);
sumV_radar = sum(V_radar,3);

SNR = zeros(1,K);
SNR_comm = zeros(1,K);
SNR_radar = zeros(1,K);

for k = 2:K + 1
    SNR(k) = h(:,:,k - 1)'*V(:,:,k)*h(:,:,k - 1)/(h(:,:,k - 1)'*(sumV - V(:,:,k))*h(:,:,k - 1) + sigma2);
    SNR_comm(k) = h(:,:,k - 1)'*V_comm(:,:,k)*h(:,:,k - 1)/(h(:,:,k - 1)'*(sumV_comm - V_comm(:,:,k))*h(:,:,k - 1) + sigma2);
    SNR_radar(k) = h(:,:,k - 1)'*V_radar(:,:,k)*h(:,:,k - 1)/(h(:,:,k - 1)'*(sumV_radar - V_radar(:,:,k))*h(:,:,k - 1) + sigma2);
end

    %% beampattern gain %%
p = zeros(1,Ntheta); % beampattern gains
aux = zeros(N,N);

for k = 1:K + 1
    aux = aux + V(:,:,k);
end

for i = 1:Ntheta
    p(i) = a(:,i)'*aux*a(:,i);
end

%% figures %%
figure();
plot(theta*180/pi,abs(p),'b','LineWidth',1);
hold on;
grid on;
plot([theta(1)*180/pi theta(end)*180/pi],[gamma(end) gamma(end)],':k','LineWidth',1);
plot(sensing_directions*180/pi,zeros(1,M),'xr','MarkerSize',10,'LineWidth',2);
plot(phi*180/pi,zeros(1,K),'o','Color',1/255*[0,104,87],'MarkerSize',10,'LineWidth',2);
for m = 1:M
    plot([sensing_directions(m)*180/pi sensing_directions(m)*180/pi],[0 1],'-.r','LineWidth',1);
end
for k = 1:K
    plot([phi(k)*180/pi phi(k)*180/pi],[0 1],'-.','Color',1/255*[0,104,87],'LineWidth',1);
end
hold off;
xlim([-90 90]);
ylim([0 0.2]);
legend('ISAC signal (N = 16)','\Gamma threshold','Target directions','User directions');
xlabel('\theta (degrees)');
ylabel('Beampattern gain');

figure();
plot(10*log10(gamma/0.001),EE_opt,'-or','LineWidth',1);
grid on;
hold on;
plot(10*log10(gamma/0.001),EE_comm*ones(1,6),'-ob','LineWidth',1); 
plot(10*log10(gamma/0.001),EE_radar*ones(1,6),'-o','Color',1/255*[0,104,87],'LineWidth',1); 
xlabel('\Gamma (dBm)');
ylabel('\eta (bps/J/Hz)');
legend('Proposed design','Communication only design','Sensing dominated design');
ylim([0 12]);
