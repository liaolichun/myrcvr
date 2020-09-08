function [cohRes] = AcqusitionByFFT(inpDat, fc, svnum, mscnt,start,GB)
%% 利用FFT方法实现伪码相位纬度的并行搜索
%% inpDat, 输入数据（已经降采样到2.048MHz或4.096MHz，前者为GPS，后者为BD）
%% fc, 降采样后的中频频率
%% svnum,伪随机码号
%% mscnt, 相干长度（以毫秒为单位）
%% start，起始点数
%% GB, 输入数据的格式， 1为BD，0为GPS

global gpsCode2048;
global bdCode4096;

%% 设置多普勒搜索范围
freq_low = -5000;
freq_high = 5000;

%% 重采样本地伪随机码
if GB==1
    fs_rate = 4.096E6 ;% code chip sampling rate
    nn = mscnt*4096;   % points of code samples
    code = repmat(bdCode4096(svnum,:),1,mscnt);  %for BD , only do 2ms coherent int
    freq_step = 1000/mscnt;  %freq resolution when doing closer lookup after code wide-off
    if svnum<6  %GEO, doppler range shrinks
        freq_low = -2000;
        freq_high = 2000;  
    end;
    dopp_bin = [freq_low:freq_step:freq_high];
else
    fs_rate = 2.048E6; 
    nn = mscnt*2048;
    code = repmat(gpsCode2048(svnum,:),1,mscnt);
%     code = repmat(gpsCode2048(svnum,:),1,mscnt*2);%llc
    freq_step = 1000/mscnt;  %freq resolution when doing closer lookup after code wide-off
    dopp_bin = [freq_low:freq_step:freq_high];
end;

%% FFT进行信号捕获过程
codeRange = [0:nn-1];
cohRes = zeros(length(dopp_bin), nn);
ts = 1/fs_rate;
xf = fft(inpDat(start:start-1+nn));
% start=500;%llc
%  xf = fft(code((start:start-1+nn)));%llc
fr = fc + freq_high ;
lc = code.*exp(j*2*pi*fr*ts*codeRange);
lcf = fft(lc);
for i = 1:length(dopp_bin),
    tmpclf = circshift(lcf',-i);
    cohRes(i,:) = ifft(xf.*conj(tmpclf'));
end






