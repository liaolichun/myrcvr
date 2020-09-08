clc;clear;close all;

% global gpsCode2048;
global_init();
IfFreq = 2e6;
Fs = 19.2e6;
data_per_byte = 4;

Ms_cnt = 4;%unit ms
BufferSize = ceil((Ms_cnt*1e-3)*Fs*(1/data_per_byte));%each byte 4 data unit byte
file = fopen('./if2e6fs19p2e6.dat','r');%MSB->LSB:Gdata1->Gdata4
[SigBuf,count] = fread(file,BufferSize,'uint8');

IfData = zeros(1,data_per_byte*count);
for n=1:count,
   IfData(4*n-3) = Convertbit(SigBuf(n),3);%sigbuf[3 2 1 0]->data[1 2 3 4]
   IfData(4*n-2) = Convertbit(SigBuf(n),2);
   IfData(4*n-1) = Convertbit(SigBuf(n),1);
   IfData(4*n) = Convertbit(SigBuf(n),0);
end

DownSampleRate = 2.048e6;
[I_Data,Q_Data] = Downsample(IfData,size(IfData,2),DownSampleRate,0,Fs,IfFreq);
Ms_cnt = 1;
for sv=2,
%% 基于FFT的方法,doppler和cdph和时域的有差异
% gpsCode2048 = zeros(33,2048);
% for i=1:33,
%     gpsPrnTbl(i,:) = GpsCodeGen(i);
% end;
% for i=1:2048,
%     k=floor(i*1.023E6/2.048E6);
%     for m=1:33
%         gpsCode2048(m,i) = gpsPrnTbl(m,( mod(k,1023)+1 ) );
%     end;
% end;
% gpsCoherentRes = AcquisitionByFFT(complex(I_Data,Q_Data), 0, sv, Ms_cnt,1,0);
% yym= abs(gpsCoherentRes(:,1:2048));
% [dop_n,n]= size(yym);
% [amp crw] = max(max(yym'));
% [amp crn] = max(max(yym));
% magMean = mean(mean(yym));
% figure();
% [x, y]=meshgrid(1:n, 1:dop_n);
% mesh(x,y,yym);
% title(['FFT svnum= ',num2str(sv)]);
% acqres.doppler = 5000 - (crw-1)*1000/Ms_cnt;
% acqres.cdph = 2048+2 - crn
%% 基于时域的方法，特别慢            
    CAcode1023 = cacode_gen(sv);
    resample_CAcode = resample(CAcode1023,DownSampleRate,Ms_cnt);

%     StartFreq = -5000;
%     FreqStep = 1000/Ms_cnt;
%     EndFreq = -4000;   
%     res_pool = DoTimeDominateAcq(I_Data,Q_Data,resample_CAcode,StartFreq,...
%         EndFreq,FreqStep,DownSampleRate,Ms_cnt);
%     %%
%     [amp,cdph] = max(max(res_pool));
%     [amp,doppler] = max(max(res_pool.'));
%     magMean = mean(mean(res_pool));
%     figure();
%     [dop,nn] = size(res_pool);
%     [x, y]=meshgrid(1:nn, 1:dop);
%     mesh(x,y,res_pool);
%     title(['时域 svnum= ',num2str(sv)]);
%     %% 
%     acqres.doppler = -((doppler - 1)*FreqStep + StartFreq);%4500
%     acqres.cdph = cdph;%875
    acqres.doppler = 4500;
    acqres.cdph = 875;

    %early early
    I_InData = I_Data(2048-acqres.cdph:2048-acqres.cdph+2048*Ms_cnt);
    Q_InData = Q_Data(2048-acqres.cdph:2048-acqres.cdph+2048*Ms_cnt);
    CarrierNcoStep = -acqres.doppler/DownSampleRate * 2^32;
    if CarrierNcoStep < 0,
        CarrierNcoStep = 2^32 + CarrierNcoStep;
    end
    CarrierNcoPhase = 0;
    CodeNcoStep = 1.023e6/DownSampleRate * 2^32;
    CodeNcoPhase = 0;
    [I_E,Q_E] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhase,CodeNcoPhase,Ms_cnt);
    sqrt(I_E^2+Q_E^2)
    %early
    I_InData = I_Data(2048-acqres.cdph+1:2048-acqres.cdph+1+2048*Ms_cnt);
    Q_InData = Q_Data(2048-acqres.cdph+1:2048-acqres.cdph+1+2048*Ms_cnt);
    CarrierNcoStep = -acqres.doppler/DownSampleRate * 2^32;
    if CarrierNcoStep < 0,
        CarrierNcoStep = 2^32 + CarrierNcoStep;
    end
    CarrierNcoPhase = 0;
    CodeNcoStep = 1.023e6/DownSampleRate * 2^32;
    CodeNcoPhase = 0;
    [I_E,Q_E] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhase,CodeNcoPhase,Ms_cnt);
    sqrt(I_E^2+Q_E^2)
    %prompt I_Data eq.   2048-start = (acqres.cdph-1) - 1
    I_InData = I_Data(2048-acqres.cdph+2:2048-acqres.cdph+2+2048*Ms_cnt);
    Q_InData = Q_Data(2048-acqres.cdph+2:2048-acqres.cdph+2+2048*Ms_cnt);
    CarrierNcoStep = -acqres.doppler/DownSampleRate * 2^32;
    if CarrierNcoStep < 0,
        CarrierNcoStep = 2^32 + CarrierNcoStep;
    end
    CarrierNcoPhase = 0;
    CodeNcoStep = 1.023e6/DownSampleRate * 2^32;
    CodeNcoPhase = 0;
    [I_P,Q_P] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhase,CodeNcoPhase,Ms_cnt);
    sqrt(I_P^2+Q_P^2)
    
    
    
    %late
    I_InData = I_Data(2048-acqres.cdph+3:2048-acqres.cdph+3+2048*Ms_cnt);
    Q_InData = Q_Data(2048-acqres.cdph+3:2048-acqres.cdph+3+2048*Ms_cnt);
    CarrierNcoStep = -acqres.doppler/DownSampleRate * 2^32;
    if CarrierNcoStep < 0,
        CarrierNcoStep = 2^32 + CarrierNcoStep;
    end    
    CarrierNcoPhase = 0;
    CodeNcoStep = 1.023e6/DownSampleRate * 2^32;
    CodeNcoPhase = 0;
    [I_L,Q_L] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhase,CodeNcoPhase,Ms_cnt);
    sqrt(I_L^2+Q_L^2)
    
    %late late
    I_InData = I_Data(2048-acqres.cdph+4:2048-acqres.cdph+4+2048*Ms_cnt);
    Q_InData = Q_Data(2048-acqres.cdph+4:2048-acqres.cdph+4+2048*Ms_cnt);
    CarrierNcoStep = -acqres.doppler/DownSampleRate * 2^32;
    if CarrierNcoStep < 0,
        CarrierNcoStep = 2^32 + CarrierNcoStep;
    end    
    CarrierNcoPhase = 0;
    CodeNcoStep = 1.023e6/DownSampleRate * 2^32;
    CodeNcoPhase = 0;
    [I_L,Q_L] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhase,CodeNcoPhase,Ms_cnt);
    sqrt(I_L^2+Q_L^2)
end