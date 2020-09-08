function ifdata = Convertbit( SigBuff,index )
%CONVERT 此处显示有关此函数的摘要
%   此处显示详细说明
    bitmap = [1 3 -1 -3];
    tmp = bitand(bitshift(SigBuff,-2*index),3);
    ifdata = bitmap(tmp+1);
end

