clear all;
close all;

N = 16; % number of antennas at the BS
K = 2; % number of users
M = 4; % number of sensing directions

s0 = zeros(N,1); % dedicated radar signal
s = round(randn(K,1)); % data signals
V = zeros(N,K); % matrix containing the beamforming vectors for all users
x = V*s + s0; % transmited signal

V0 = toeplitz(xcorr(s0)); % dedicated radar signal autocorrelation matrix

h = zeros(N,K); % matrix containing the channels from the BS to all users




