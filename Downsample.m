function [ I_Data Q_Data ] = Downsample( inputdata,length,samplerate,fc,iffs,iffreq )
%DOWNSAMPLE 此处显示有关此函数的摘要
%   此处显示详细说明
global sin_table;
global cos_table;
    cnt = 1;
    carrierNCOphase = 0;
    carrierStep = round(iffreq / iffs * 2^32);
    samplerateStep = round(samplerate / iffs * 2^32);
    sampleNCOphase = 0;
    sampleNCOphase_old = 0;
    I_Data(1) = 0;
    Q_Data(1) = 0;
    for n=1:1:length,        
        carrierNCOphase = mod(carrierNCOphase + carrierStep,2^32-1);
        phase_idx = bitand(bitshift(carrierNCOphase,-28),15)+1;
        I_Data(cnt) = I_Data(cnt) + inputdata(n)*cos_table(phase_idx);
        Q_Data(cnt) = Q_Data(cnt) + inputdata(n)*-sin_table(phase_idx);
        
        sampleNCOphase = mod(sampleNCOphase + samplerateStep,2^32-1);
        if(sampleNCOphase < sampleNCOphase_old)
            cnt = cnt + 1;
            I_Data(cnt) = 0;
            Q_Data(cnt) = 0;
        end
        sampleNCOphase_old = sampleNCOphase;
    end
end

