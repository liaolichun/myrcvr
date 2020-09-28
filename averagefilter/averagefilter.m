clc;clear;close all;
l = [2 3 4 5 6 7];
% l = 3;
N = length (l);

for i = 1:N
    L = l (i);
    B = (1/L)* ones (1, L);
    [H, F] = freqz (B, 1, 100);
    A (i,:) = abs (H);
end;
plot(A(1,:),'r');
hold on;
plot(A(2,:),'g');
plot(A(3,:),'b');
plot(A(4,:),'y');
plot(A(5,:),'m');
plot(A(6,:),'c');