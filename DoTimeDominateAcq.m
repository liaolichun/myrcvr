function [res_pool] = DoTimeDominateAcq( I_Data,Q_Data,resample_CAcode,...
    startbin,endbin,bin_step,fs,ms_cnt)
%DOTIMEDOMINATEACQ 此处显示有关此函数的摘要
%   此处显示详细说明
global sin_table;
global cos_table;

carrierNcoPhase = 0;
res_pool = zeros((endbin-startbin)/bin_step + 1,2048);
  
for freq=startbin:bin_step:endbin,
    carrierNcoStep = freq/fs * 2^32;
    if carrierNcoStep < 0,
        carrierNcoStep = 2^32 + carrierNcoStep;
    end
    for hyp=1:1:2048,     
        tmp_CAcode = [resample_CAcode(hyp:end) resample_CAcode(1:hyp-1)]; 
        integrate_I = 0;
        integrate_Q = 0;
        for nn=1:1:2048*ms_cnt,            
            carrierNcoPhase = mod(carrierNcoPhase + carrierNcoStep,2^32-1);
            phase_idx = bitand(bitshift(carrierNcoPhase,-28),15)+1;
            buf_I = I_Data(nn) * cos_table(phase_idx) - Q_Data(nn) * sin_table(phase_idx);
            buf_Q = I_Data(nn) * sin_table(phase_idx) + Q_Data(nn) * cos_table(phase_idx);
            integrate_I = integrate_I + buf_I * tmp_CAcode(nn);
            integrate_Q = integrate_Q + buf_Q * tmp_CAcode(nn);            
        end
        res_pool((freq-startbin)/bin_step+1,hyp) = sqrt(integrate_I^2+integrate_Q^2);
        carrierNcoPhase = 0;
    end
end
% [amp,cdph] = max(max(res_pool));
% [amp,doppler] = max(max(res_pool.'));
% magMean = mean(mean(res_pool));
% figure();
% [dop,nn] = size(res_pool);
% [x, y]=meshgrid(1:nn, 1:dop);
% mesh(x,y,res_pool);
end

