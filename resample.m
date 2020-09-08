function [ res_CAcode ] = resample( CAcode,DownSampleRate,ms_cnt );
%RESAMPLE 此处显示有关此函数的摘要
%   此处显示详细说明

    ca_freq = 1.023e6;
    SampleRateNcoStep = ca_freq/DownSampleRate * 2^32;
    SampleRateNcoPhase = 0;
    SampleRateNcoPhase_old = 0;
    cnt_res = 1;
    cnt_orig = 1;
    for n =1:1:DownSampleRate*(1e-3)
        res_CAcode(cnt_res) = CAcode(cnt_orig);
        cnt_res = cnt_res + 1;
        SampleRateNcoPhase = mod(SampleRateNcoPhase + SampleRateNcoStep,2^32);
        if(SampleRateNcoPhase<SampleRateNcoPhase_old)
            cnt_orig = cnt_orig + 1;
        end
        SampleRateNcoPhase_old = SampleRateNcoPhase;
    end
    res_CAcode = repmat(res_CAcode,1,ms_cnt);
end

