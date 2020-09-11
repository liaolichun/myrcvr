clc;clear;close all;

global gpsCode2048;
global LUT;
global_init();
% IfFreq = 2e6;
% Fs = 19.2e6;

file = fopen('./if3p996e6fs16p369e6.dat','r');%MSB->LSB:Gdata1->Gdata4
IfFreq = 3.996E6;
Fs = 16.369E6;
data_per_byte = 2;
read_Ms_cnt = 500;%unit ms
BufferSize = ceil((read_Ms_cnt*1e-3)*Fs*(1/data_per_byte));%each byte 4 data unit byte
[SigBuf,count] = fread(file,BufferSize,'uint8');

IfData = zeros(1,data_per_byte*count);
for n=1:count,
%    IfData(4*n-3) = Convertbit(SigBuf(n),3);%sigbuf[3 2 1 0]->data[1 2 3 4]
%    IfData(4*n-2) = Convertbit(SigBuf(n),2);
%    IfData(4*n-1) = Convertbit(SigBuf(n),1);
%    IfData(4*n) = Convertbit(SigBuf(n),0);
%     IfData(2*n-1) = Convertbit(SigBuf(n),2);
%     IfData(2*n) = Convertbit(SigBuf(n),0);
    IfData(2*n-1) = LUT(SigBuf(n)+1,2);
    IfData(2*n) = LUT(SigBuf(n)+1,4);
end

DownSampleRate = 2.048e6;
[I_Data,Q_Data] = Downsample(IfData,size(IfData,2),DownSampleRate,0,Fs,IfFreq);

for sv=5,
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
% Ms_cnt = 1;
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
%     CAcode1023 = cacode_gen(sv);
%     Ms_cnt = 2;
%     resample_CAcode = resample(CAcode1023,DownSampleRate,Ms_cnt);
% 
%     StartFreq = -5000;
%     FreqStep = 1000/Ms_cnt;
%     EndFreq = 5000;   
%     res_pool = DoTimeDominateAcq(I_Data,Q_Data,resample_CAcode,StartFreq,...
%         EndFreq,FreqStep,DownSampleRate,Ms_cnt);
% %
%     [amp,cdph] = max(max(res_pool));
%     [amp,doppler] = max(max(res_pool.'));
%     magMean = mean(mean(res_pool));
%     figure();
%     [dop,nn] = size(res_pool);
%     [x, y]=meshgrid(1:nn, 1:dop);
%     mesh(x,y,res_pool);
%     title(['时域 svnum= ',num2str(sv)]);
% % 
%     acqres.doppler = -((doppler - 1)*FreqStep + StartFreq);%-2500
%     acqres.cdph = cdph;%817
%% Tracking base the result of timedomain 
%% initial sv2 -2500 817 sv5 -1500 1204  PULLIN
    acqres.doppler = -1500;%-2500
    acqres.cdph = 1204;%817;

    total_ms = read_Ms_cnt - 1;
    freq_change = 0;
    cdph_change = 0;
    
    CAcode1023 = cacode_gen(sv);
    CarrierNcoStep = -acqres.doppler/DownSampleRate * 2^32;
    if CarrierNcoStep < 0,
        CarrierNcoStep = 2^32 + CarrierNcoStep;
    end    
    CodeNcoStep = (1.023e6+acqres.doppler/1540)/DownSampleRate * 2^32;
    
    CarrierNcoPhaseE = 0;
    CarrierNcoPhaseP = 0;
    CarrierNcoPhaseL = 0;
    CodeNcoPhaseE = 0;
    CodeNcoPhaseP = 0;
    CodeNcoPhaseL = 0;
    codeidxE = 0;
    codeidxP = 0;
    codeidxL = 0;
    delta_cdph_old = 0;
    Ip_old = 0;
    Qp_old = 0;
    ms_counter = 1;% range 1-20
    bitedge = zeros(1,20);
    delta_freq_old = 0;
    cnta = 1;
    cntb = 1;
    cntc = 1;
    cntd = 1;
    for i=1:1:total_ms
        CarrierNcoStep = CarrierNcoStep-freq_change/DownSampleRate * 2^32;
        CodeNcoStep = CodeNcoStep + cdph_change/DownSampleRate * 2^32;
        Ms_cnt = 1;
%% early
        I_InData = I_Data(2048-acqres.cdph+1:2048-acqres.cdph+1+2048*Ms_cnt);
        Q_InData = Q_Data(2048-acqres.cdph+1:2048-acqres.cdph+1+2048*Ms_cnt);
        [Ie,Qe,codeNcoE,carrierNcoE,CodeIdxE] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhaseE,CodeNcoPhaseE,Ms_cnt,0,0,codeidxE);
        codeidxE = CodeIdxE;
        CarrierNcoPhaseE = carrierNcoE;
        CodeNcoPhaseE = codeNcoE;
%% prompt I_Data eq.   2048-start = (acqres.cdph-1) - 1        
        I_InData = I_Data(2048-acqres.cdph+2:2048-acqres.cdph+2+2048*Ms_cnt);
        Q_InData = Q_Data(2048-acqres.cdph+2:2048-acqres.cdph+2+2048*Ms_cnt);
        [Ip,Qp,codeNcoP,carrierNcoP,CodeIdxP] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhaseP,CodeNcoPhaseP,Ms_cnt,0,0,codeidxP);
        codeidxP = CodeIdxP;
        CarrierNcoPhaseP = carrierNcoP;
        CodeNcoPhaseP = codeNcoP;
%% late
        I_InData = I_Data(2048-acqres.cdph+3:2048-acqres.cdph+3+2048*Ms_cnt);
        Q_InData = Q_Data(2048-acqres.cdph+3:2048-acqres.cdph+3+2048*Ms_cnt);
        [Il,Ql,codeNcoL,carrierNcoL,CodeIdxL] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhaseL,CodeNcoPhaseL,Ms_cnt,0,0,codeidxL);
        codeidxL = CodeIdxL;
        CarrierNcoPhaseL = carrierNcoL;
        CodeNcoPhaseL = codeNcoL;
%% code phase discriminator
        SE = sqrt(Ie^2+Qe^2);
        SP = sqrt(Ip^2+Qp^2);
        SL = sqrt(Il^2+Ql^2);
        delta_cdph(i) = (SE - SL)/SP;        
        cdph_change = delta_cdph(i)*0.5 + (delta_cdph(i) - delta_cdph_old)*2;
        delta_cdph_old = delta_cdph(i);
%% freq discriminator
        dot_p = Ip_old * Ip + Qp_old * Qp;
        cross_p = Ip_old * Qp - Qp_old * Ip;
        delta_freq(i) = atan2(cross_p,dot_p)/(2*pi*Ms_cnt*1e-3); 
        if abs(delta_freq(i) - delta_freq_old)<250
            freq_change = delta_freq(i) * 0.05;        
        end
        delta_freq_old = delta_freq(i);
        theta_1(cnta) = atan2(Qp,Ip);
        cnta = cnta+1;
        Ip_old = Ip;
        Qp_old = Qp;
%%
        I_Data = I_Data(2048+1:end);
        Q_Data = Q_Data(2048+1:end);
        ms_counter = mod(ms_counter,20) + 1;
    end
%     figure(1);
%     plot(delta_cdph);
%     figure(2);
%     plot(delta_freq);
%     figure(3);
%     plot(theta_1);
%% get new data for TRACKING
    phd_old = 0;
    bitedge_done = 0;
    I_buffer = zeros(1,total_ms);
    Q_buffer = zeros(1,total_ms);
    BN = 0.5;
    c1 = (BN/0.53)^2;
    c2 = 1.414*BN/0.53;
    while 1
        read_Ms_cnt = 1000;%unit ms
        BufferSize = ceil((read_Ms_cnt*1e-3)*Fs*(1/data_per_byte));%each byte 4 data unit byte
        [SigBuf,count] = fread(file,BufferSize,'uint8');
        if count < BufferSize
            break;
        end
        IfData = zeros(1,data_per_byte*count);
        for n=1:count,
%             IfData(2*n-1) = Convertbit(SigBuf(n),2);
%             IfData(2*n) = Convertbit(SigBuf(n),0);
            IfData(2*n-1) = LUT(SigBuf(n)+1,2);
            IfData(2*n) = LUT(SigBuf(n)+1,4);
        end
        [I_tmp,Q_tmp] = Downsample(IfData,size(IfData,2),DownSampleRate,0,Fs,IfFreq);
        I_Data = [I_Data I_tmp];
        Q_Data = [Q_Data Q_tmp];
        total_ms = read_Ms_cnt;
        for i=1:1:total_ms
            CarrierNcoStep = CarrierNcoStep-freq_change/DownSampleRate * 2^32;
            CodeNcoStep = CodeNcoStep + cdph_change/DownSampleRate * 2^32 + (freq_change/1540)/DownSampleRate * 2^32;
%% early
            I_InData = I_Data(2048-acqres.cdph+1:2048-acqres.cdph+1+2048*Ms_cnt);
            Q_InData = Q_Data(2048-acqres.cdph+1:2048-acqres.cdph+1+2048*Ms_cnt);
            [Ie,Qe,codeNcoE,carrierNcoE,CodeIdxE] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhaseE,CodeNcoPhaseE,Ms_cnt,0,0,codeidxE);
            codeidxE = CodeIdxE;
            CarrierNcoPhaseE = carrierNcoE;
            CodeNcoPhaseE = codeNcoE;
%% prompt I_Data eq.   2048-start = (acqres.cdph-1) - 1        
            I_InData = I_Data(2048-acqres.cdph+2:2048-acqres.cdph+2+2048*Ms_cnt);
            Q_InData = Q_Data(2048-acqres.cdph+2:2048-acqres.cdph+2+2048*Ms_cnt);
            [Ip,Qp,codeNcoP,carrierNcoP,CodeIdxP] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhaseP,CodeNcoPhaseP,Ms_cnt,0,0,codeidxP);
            codeidxP = CodeIdxP;
            CarrierNcoPhaseP = carrierNcoP;
            CodeNcoPhaseP = codeNcoP;
%% late
            I_InData = I_Data(2048-acqres.cdph+3:2048-acqres.cdph+3+2048*Ms_cnt);
            Q_InData = Q_Data(2048-acqres.cdph+3:2048-acqres.cdph+3+2048*Ms_cnt);
            [Il,Ql,codeNcoL,carrierNcoL,CodeIdxL] = CorherentSum(I_InData,Q_InData,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhaseL,CodeNcoPhaseL,Ms_cnt,0,0,codeidxL);
            codeidxL = CodeIdxL;
            CarrierNcoPhaseL = carrierNcoL;
            CodeNcoPhaseL = codeNcoL;
%% code phase discriminator
            
%% cdph discriminator        
            if bitedge_done == 0 %no bitsync 
                SE = sqrt(Ie^2+Qe^2);% code loop
                SP = sqrt(Ip^2+Qp^2);
                SL = sqrt(Il^2+Ql^2);
                delta_cdph(i) = (SE - SL)/SP;        
                cdph_change = delta_cdph(i)*0.02 + (delta_cdph(i) - delta_cdph_old)*0.1;
                delta_cdph_old = delta_cdph(i);
                
                dot_p = Ip_old * Ip + Qp_old * Qp;%carrier loop
                cross_p = Ip_old * Qp - Qp_old * Ip;
                delta_freq(i) = atan2(cross_p,dot_p)/(2*pi*Ms_cnt*1e-3);
                if abs(delta_freq(i) - delta_freq_old)<250
                    freq_change = delta_freq(i) * 0.005 ;
                end
                delta_freq_old = delta_freq(i);
                theta_2(cntb) = atan2(Qp,Ip);
                cntb = cntb + 1;
%                 figure(4);
%                 plot(delta_cdph);
            else                
                sum_20ms_Ie = sum_20ms_Ie + Ie;
                sum_20ms_Qe = sum_20ms_Qe + Qe;
                sum_20ms_Ip = sum_20ms_Ip + Ip;
                sum_20ms_Qp = sum_20ms_Qp + Qp;
                sum_20ms_Il = sum_20ms_Il + Il;
                sum_20ms_Ql = sum_20ms_Ql + Ql;
                sum_10ms_Ie = sum_10ms_Ie + Ie;
                sum_10ms_Qe = sum_10ms_Qe + Qe;
                sum_10ms_Ip = sum_10ms_Ip + Ip;
                sum_10ms_Qp = sum_10ms_Qp + Qp;
                sum_10ms_Il = sum_10ms_Il + Il;
                sum_10ms_Ql = sum_10ms_Ql + Ql;
                if mod(ms_counter,10) == 0   %10ms carrier loop
                    if sum_10ms_Ip > 0
                        phd(t_cnt) = sum_10ms_Qp/sqrt(sum_10ms_Ip^2+sum_10ms_Qp^2);
                    else
                        phd(t_cnt) = -sum_10ms_Qp/sqrt(sum_10ms_Ip^2+sum_10ms_Qp^2);
                    end
%                     freq_change = (phd(t_cnt) - phd_old)*0.04+ ...
%                         phd(t_cnt) * 0.002;
                    freq_change = (phd(t_cnt) - phd_old)*c2+ ...
                        phd(t_cnt) * c1 * 0.01;
                    phd_old = phd(t_cnt);                    
                    t_cnt = t_cnt + 1;
                    sum_10ms_Ie = 0;
                    sum_10ms_Qe = 0;
                    sum_10ms_Ip = 0;
                    sum_10ms_Qp = 0;
                    sum_10ms_Il = 0;
                    sum_10ms_Ql = 0;
                    freq_chage_tmp (cntd) = freq_change;
                    cntd = cntd + 1;
                    figure(5);
                    plot(phd);
                    title('phd');
                    figure(6);
                    plot(freq_chage_tmp);
                    title('freq chage tmp');
                end
                if mod(ms_counter,20) == 0  %20ms code loop
                    magE20ms = sqrt(sum_20ms_Ie^2+sum_20ms_Qe^2);
                    magL20ms = sqrt(sum_20ms_Il^2+sum_20ms_Ql^2);
                    delta_cdph2(cntc) = 0.5* (magE20ms - magL20ms) / (magE20ms + magL20ms);
                    cdph_change = delta_cdph2(cntc)*0.2 + (delta_cdph2(cntc) - delta_cdph_old)*0.8;    
%                     cdph_change = delta_cdph2(cntc)*0.2;
                    delta_cdph_old = delta_cdph2(cntc);
                    cntc = cntc + 1;
                    sum_20ms_Ie = 0;
                    sum_20ms_Qe = 0;
                    sum_20ms_Ip = 0;
                    sum_20ms_Qp = 0;
                    sum_20ms_Il = 0;
                    sum_20ms_Ql = 0;
                    figure(7);
                    plot(delta_cdph2);
                    title('delta cdph2');
                end
                I_buffer = [I_buffer(2:end) Ip];
                Q_buffer = [Q_buffer(2:end) Qp];
                figure(8);
                plot(I_buffer,'r');
                hold on;
                plot(Q_buffer,'b');
                title('I red,Q blue');
                hold off;
                pause(0.001);
            end
            Ip_old = Ip;
            Qp_old = Qp;
%% push data in and do bitsync
            I_Data = I_Data(2048*Ms_cnt+1:end);
            Q_Data = Q_Data(2048*Ms_cnt+1:end);
            if dot_p < 0 && i<501 && bitedge_done == 0
                bitedge(ms_counter) = bitedge(ms_counter) + 1;
            end
            if i == 501 && bitedge_done == 0
                [amp mscounter] = max(bitedge);%TODO
                ms_counter = ms_counter - (mscounter - 1);
                bitedge_done = 1;
                sum_20ms_Ie = 0;
                sum_20ms_Qe = 0;
                sum_20ms_Ip = 0;
                sum_20ms_Qp = 0;
                sum_20ms_Il = 0;
                sum_20ms_Ql = 0;
                sum_10ms_Ie = 0;
                sum_10ms_Qe = 0;
                sum_10ms_Ip = 0;
                sum_10ms_Qp = 0;
                sum_10ms_Il = 0;
                sum_10ms_Ql = 0;
                t_cnt = 1;
            end  
            ms_counter = mod(ms_counter,20) + 1;
        end
%         figure(1);
%         plot(delta_cdph(1:total_ms));
%         figure(2);
%         plot(delta_freq);
%         figure(3);
%         plot(IQ_theta);
%         pause(0.0001);
    end
end