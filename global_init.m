function []=global_init();
global sin_table;
global cos_table;
global LUT;

sin_table = [0,6,11,15,16,15,11, 6, 0,-6,-11,-15,-16,-15,-11,-6];
cos_table = [16,15,11, 6, 0,-6,-11,-15,-16,-15,-11,-6, 0, 6,11,15];
valuemap = [1 3 -1 -3];
LUT = zeros(256,4);

for i = 0:1:255
   LUT(i+1,1) = valuemap(bitand(bitshift(i,-6),3)+1);
   LUT(i+1,2) = valuemap(bitand(bitshift(i,-4),3)+1);
   LUT(i+1,3) = valuemap(bitand(bitshift(i,-2),3)+1);
   LUT(i+1,4) = valuemap(bitand(bitshift(i,0),3)+1);
end
    
end