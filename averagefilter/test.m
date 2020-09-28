clc;clear;close all;
% w = 0:0.01:2*pi;
% H = 1+ exp(-j*w) + exp(-j*2*w);
% plot(abs(H));
% 
f = 0:0.01:1;
H = 1 + exp(-j*2*pi*f) + exp(-j*2*2*pi*f) ;%+ exp(-j*3*2*pi*f) + exp(-j*4*2*pi*f)...
   % + exp(-j*5*2*pi*f) + exp(-j*6*2*pi*f) + exp(-j*7*2*pi*f);
plot(abs(H));

f1 = 30;
f2 = 200;
fs = 600;
t = 1:1:300;
s = exp(j*2*pi*f1*t/fs) + exp(j*2*pi*f2*t/fs);
F_s = fft(s);
figure(1);
plot(abs(F_s));

filter_s = s(3:end) + s(2:end-1) + s(1:end-2);
F_filter_s = fft(filter_s);
figure(2);
plot(abs(F_filter_s));