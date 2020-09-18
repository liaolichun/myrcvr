clc;clear;close all;

f1 = 1500;
f2 = 1490;
fs = 2.048e6;
phi = pi/4;

time_ms = 1000;

n = time_ms*1e-3*fs;
t = 1:1:n;
SNR = 30;% unit db
noise_power = 10^-(SNR/10);
n_I = 0.707*noise_power*randn(1,n);
n_Q = 0.707*noise_power*randn(1,n);
ifdata = exp(j*(2*pi*f1/fs*t+phi))+complex(n_I,n_Q);phd_old = 0;
BB_data = zeros(1,time_ms);
BN = 16;
c1 = (BN/0.53)^2;
c2 = 1.414*BN/0.53;
local_phase = 0;

integrate_time = 5; %unit ms
buffer_I = zeros(1,time_ms/integrate_time);
buffer_Q = zeros(1,time_ms/integrate_time);
for m =1:1:time_ms/integrate_time
    for x =1:1:1e-3*fs*integrate_time
        local_phase = local_phase + 2*pi*f2/fs;
        BB_data(m) = BB_data(m) + ifdata((m-1)*integrate_time*fs*1e-3 + x) * exp(-j*local_phase);
    end    
    buffer_I =  [buffer_I(2:end) real(BB_data(m))];
    buffer_Q =  [buffer_Q(2:end) imag(BB_data(m))];
%     phd(m) = angle(BB_data(m));  
    if real(BB_data(m)) > 0
        phd(m) = imag(BB_data(m))/sqrt(real(BB_data(m))^2+imag(BB_data(m))^2);
    else
        phd(m) = -imag(BB_data(m))/sqrt(real(BB_data(m))^2+imag(BB_data(m))^2);
    end
    freq_change = phd(m)*integrate_time*1e-3*c1 + (phd(m) - phd_old)*c2;%
    f2 = f2 + freq_change;
    phd_old = phd(m);
    figure(1);
    plot(buffer_I,'r');
    hold on;
    plot(buffer_Q,'b');
    hold off;
    title('I red,Q blue');
    figure(2);
    plot(phd);
    pause(0.0001);
end
