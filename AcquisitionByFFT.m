function [cohRes] = AcqusitionByFFT(inpDat, fc, svnum, mscnt,start,GB)
%% ����FFT����ʵ��α����λγ�ȵĲ�������
%% inpDat, �������ݣ��Ѿ���������2.048MHz��4.096MHz��ǰ��ΪGPS������ΪBD��
%% fc, �����������ƵƵ��
%% svnum,α������
%% mscnt, ��ɳ��ȣ��Ժ���Ϊ��λ��
%% start����ʼ����
%% GB, �������ݵĸ�ʽ�� 1ΪBD��0ΪGPS

global gpsCode2048;
global bdCode4096;

%% ���ö�����������Χ
freq_low = -5000;
freq_high = 5000;

%% �ز�������α�����
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

%% FFT�����źŲ������
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






