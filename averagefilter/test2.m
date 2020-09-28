clc;clear;close all;
fs = 500;
n = 1:1000;
s1 = sin (2*pi*10*n*(1/fs));
s2 = sin (2*pi*15*n*(1/fs));
s3 = sin (2*pi*20*n*(1/fs));
s4 = sin (2*pi*25*n*(1/fs));
s5 = sin (2*pi*30*n*(1/fs));
x = s1 + s2 + s3 + s4 + s5;
v = 2*sin (2*pi*0.2*n*(1/fs));
s = x + v;


L = 128;
b = (1/L)*ones (1, L);
y = filter (b, 1, s);
plot(y);