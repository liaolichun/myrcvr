function [ca_used] = GpsCodeGen(svnum);
%% ²úÉúBDÎ±Ëæ»úÂë£¬
%% svnum:1--32 

% gs2= [5; 6; 7; 8; 17; 18; 139; 140; 141; 251; 252; 254; 255; 256; 257; 258; 469; 470; 471; 472; 473; 474; 509; 512; 513; 514; 515; 516; 859; 860; 861; 862];
gs2= [5; 6; 7; 8; 17; 18; 139; 140; 141; 251; 252; 254; 255; 256; 257; 258; 469; 470; 471; 472; 473; 474; 509; 512; 513; 514; 515; 516; 859; 860; 861; 862;863;950;947;948;950];
g2shift = gs2(svnum,1);
 
% G1 LFSR: x^10+x^3+1
% generate G1 code
reg = -1*ones(1,10);
for i = 1:1023,
    g1(i) = reg(10);
    save1 = reg(3) * reg(10);
    reg(1,2:10) = reg(1,1:9);
    reg(1,1) = save1;
end;

% G2j LFSR: x^10+x^9+x^8+x^6+x^3+x^2+1
% genereate G2 code
reg = -1*ones(1,10);
for i = 1:1023,
    g2(i) = reg(10);
    save2 = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
    reg(1,2:10) = reg(1,1:9);
    reg(1) = save2;
end

% shift G2 code , according to ICD200 
g2tmp(1,1:g2shift) = g2(1,1023-g2shift+1:1023);
g2tmp(1,g2shift+1:1023) = g2(1, 1:1023-g2shift);
g2 = g2tmp;

ss_ca = g1.*g2;
ca_used = -ss_ca;