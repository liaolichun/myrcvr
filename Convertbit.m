function ifdata = Convertbit( SigBuff,index )
%CONVERT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    bitmap = [1 3 -1 -3];
    tmp = bitand(bitshift(SigBuff,-2*index),3);
    ifdata = bitmap(tmp+1);
end

