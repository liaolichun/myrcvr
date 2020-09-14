clc;clear;close all;

f1 = 100;
f2 = 90;
fs = 2.048e6;
phi = pi/4;

time_ms = 1000;

n = time_ms*1e-3*fs;
t = 1:1:n;
SNR = -10;% unit db
noise_power = 10^-(SNR/10);
n_I = 0.707*noise_power*randn(1,n);
n_Q = 0.707*noise_power*randn(1,n);
ifdata = exp(j*(2*pi*f1/fs*t+phi))+complex(n_I,n_Q);

phd_old = 0;
BB_data = zeros(1,time_ms);
BN = 16;
c1 = (BN/0.53)^2;
c2 = 1.414*BN/0.53;
local_phase = 0;

integrate_time = 5; %unit ms
for m =1:1:time_ms/integrate_time
    for x =1:1:1e-3*fs*integrate_time
        local_phase = local_phase + 2*pi*f2/fs;
        BB_data(m) = BB_data(m) + ifdata((m-1)*integrate_time*fs*1e-3 + x) * exp(-j*local_phase);
    end    
    phd(m) = angle(BB_data(m));      
    freq_change = phd(m)*integrate_time*1e-3*c1 + (phd(m) - phd_old)*c2;%
    f2 = f2 + freq_change;
    phd_old = phd(m);
end
plot(phd)