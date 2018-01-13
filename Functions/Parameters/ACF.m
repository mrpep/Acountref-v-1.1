% Código de Shin-Ichi Sato adaptado por Leonardo Pepino para cómputo de
% parámetros de Ando.

% Parametros de entrada:

% int = Integration Interval - Default = 1 s
% step = Running step - Default = 0.1 s
% tmax = Maximum delay range of ACF [s] - Default = 0.2 s
% rgtau1 = Longer limit of tau1 detection [ms] - Default = 10 ms

% Salida:

%RunningFactor(:,1) = indice
%RunningFactor(:,2) = vector temporal en segundos
%RunningFactor(:,3) = Phi(0)L
%RunningFactor(:,4) = Tau1L en ms
%RunningFactor(:,5) = Phi1 L
%RunningFactor(:,6) = TauE L -> Este es el deseado
%RunningFactor(:,7) = Corr L
%RunningFactor(:,8) = IndexTauE L
%RunningFactor(:,9) = 2*width_L*1000 ¿?

%Si la señal de entrada es binaural hay más indices:

%RunningFactor(:,10) = Phi(0)R
%RunningFactor(:,11) = Tau1R en ms
%RunningFactor(:,12) = Phi1 R
%RunningFactor(:,13) = TauE R
%RunningFactor(:,14) = Corr R
%RunningFactor(:,15) = IndexTauE R
%RunningFactor(:,16) = 2*width_R*1000
%RunningFactor(:,17) = Phi0_LR
%RunningFactor(:,18) = IACC -> También usarlo
%RunningFactor(:,19) = tau
%RunningFactor(:,20) = wiacc

function [RunningFactor] = ACF(y,fs,int,step,tmax,rgtau1)

y = afilter(y,fs);
tv = (1:length(y))'./fs;

rnn = (1:floor(step*fs):floor(length(y)/floor(step*fs))*floor(step*fs))';
runstep = length(rnn)-floor(((int+tmax)/step));

for z= 1:runstep;
    trg = rnn(z);
    index = (0:1/fs:tmax-1/fs)';
    
    data0_L = y(trg:trg+floor(int*fs)-1,1);
    data0_L(:) = data0_L(:)-mean(data0_L(:));
    max_norm_L = max(abs(data0_L(:)));
    data0_L(:) = data0_L(:)./max_norm_L;
    datamax_L = 10*log10(max(data0_L.^2));
    
    data0_L = y(trg:trg+floor((int+tmax)*fs)-1,1);
    data0_L(:) = data0_L(:)-mean(data0_L);
    data1_L = cat(1,data0_L(1:floor(int*fs)), zeros(size(data0_L(floor(int*fs)+1:end))));
    
    data_L = [data1_L data0_L];
    FFT_L = fft(data_L);
    CONJ_L = conj(FFT_L(:,2));
    ACF = real(ifft(conj(fft(data0_L)).*fft(data1_L)));
    Phi0_L = 10*log10(ACF(1,1));
    
    normalizer_L = zeros(floor(tmax*fs),1);
    n = 1;
    for lag = 0,
        tmp = data0_L(1:floor(int*fs)).^2;
        Phi_tmp = sum(tmp);
        Phi_t = sum(tmp);
        normalizer_L(n,:) = Phi_tmp;
        n = n+1;
    end
    for lag = 1:floor(tmax*fs)-1;
        Phi_tmp = Phi_tmp-data0_L(lag)^2+data0_L(floor(int*fs)+lag)^2;
        normalizer_L(n,:) = (Phi_t*Phi_tmp)^0.5;
        n = n+1;
    end
    
    NACF_L = ACF(1:floor(tmax*fs),1)./normalizer_L;
    acfplot_L = [index NACF_L];
    lgacfplot_L = [index 10*log10(abs(NACF_L))];
    
    if size(y,2) == 2;
        
        data0_R = y(trg:trg+floor(int*fs)-1,2);
        datamax_R = 10*log10(max(data0_R.^2));
        
        data0_R = y(trg:trg+floor((int+tmax)*fs)-1,2);
        data0_R(:) = data0_R(:)-mean(data0_L);
        data1_R = cat(1,data0_R(1:floor(int*fs)-1), zeros(size(data0_R(floor(int*fs):end))));
        data_R = [data1_R data0_R];
        FFT_R = fft(data_R);
        CONJ_R = conj(FFT_R(:,2));
        
        data = [CONJ_L(:,1).*FFT_L(:,1) CONJ_R(:,1).*FFT_R(:,1) CONJ_L(:,1).*FFT_R(:,1) CONJ_R(:,1).*FFT_L(:,1)];
        IFFT = ifft(data);
        ACF = real(IFFT(:,:));
        
        Phi0_R = 10*log10(ACF(1,2));
        
        normalizer_R = zeros(tmax*fs,1);
        n = 1;
        for lag = 0,
            tmp = data0_R(1:int*fs).*data0_R(1+lag:int*fs+lag);
            Phi_tmp = sum(tmp);
            Phi_t = sum(tmp);
            normalizer_R(n,:) = (Phi_t.*Phi_tmp).^0.5;
            n = n+1;
        end
        for lag = 1:floor(tmax*fs)-1;
            Phi_tmp = Phi_tmp - data0_R(lag).^2 + data0_R(floor(int*fs)+lag).^2;
            normalizer_R(n,:) = (Phi_t.*Phi_tmp).^0.5;
            n = n+1;
        end
        
        NACF_R = ACF(1:floor(tmax*fs),2)./normalizer_R;
        acfplot_R = [index NACF_R];
        lgacfplot_R = [index 10*log10(abs(NACF_R))];
        
        nrm = sqrt(ACF(1,1).*ACF(1,2));
        lftplt = cat(2, (0:-1000/fs:floor((-1.0*0.02*fs+1)*1000)/fs)', real(IFFT(1:floor(1.0*0.02*fs),3))./nrm);
        rgtplt = cat(2, (1000./fs:1000./fs:round(1.0*0.02*fs-1)*1000./fs)', real(IFFT(2:round(1.0*0.02*fs),4))./nrm);
        ccfplot = cat(1, flipud(lftplt), rgtplt);
        Phi0_LR = 10*log10(nrm);
        
        [tau1_R phi1_R width_R] = acffactor1_R(tmax, fs, rgtau1, acfplot_R);
        [taue_R corr_R index_taue_R reg0_R r_taue] = acffactor2_R(tmax, fs, acfplot_R, phi1_R, lgacfplot_R);
        [IACC tau wiacc] = iacffactor(fs, ccfplot);
        
    end
    
    [tau1_L phi1_L width_L] = acffactor1_L(tmax, fs, rgtau1, acfplot_L);
    [taue_L corr_L index_taue_L reg0_L r_taue] = acffactor2_L(tmax, fs, acfplot_L, phi1_L, lgacfplot_L);
    
    if size(y,2) == 1;
        RunningFactor(z,:)=[z (z-1)*step Phi0_L tau1_L*1000 phi1_L taue_L corr_L index_taue_L 2*width_L*1000];
    else
        RunningFactor(z,:)=[z (z-1)*step Phi0_L tau1_L*1000 phi1_L taue_L corr_L index_taue_L 2*width_L*1000 Phi0_R tau1_R*1000 phi1_R taue_R corr_R index_taue_R 2*width_R*1000 Phi0_LR IACC tau wiacc];
    end
    
    RunningFactor1(z,:,:) = [r_taue];
       
end

if size(y,2) == 1;
    RunningFactor(z,:)=[z (z-1)*step Phi0_L tau1_L*1000 phi1_L taue_L corr_L index_taue_L 2*width_L*1000];
else
    RunningFactor(z,:)=[z (z-1)*step Phi0_L tau1_L*1000 phi1_L taue_L corr_L index_taue_L 2*width_L*1000 Phi0_R tau1_R*1000 phi1_R taue_R corr_R index_taue_R 2*width_R*1000 Phi0_LR IACC tau wiacc];
end

RunningFactor1(z,:,:) = [r_taue];

if min(RunningFactor(:,6))<10;
    lowerrange = 1;
else lowerrange = 10;
end
%% external functions

% --------------------------------------------------------------------
function [tau1_L phi1_L width_L] = acffactor1_L(tmax, fs, rgtau1, acfplot_L)

% calculation of tau_1 and phi_1
for i = 1:floor(rgtau1/1000*fs);
    if acfplot_L(i,2) < 0
        tg = i;,break,
    end
end

for i = 1:floor(rgtau1/1000*fs);
    tc(i,:) = exp(-acfplot_L(i,1)/0.01);
end
acfplot2_L = acfplot_L(1:floor(rgtau1/1000*fs),2).*tc(1:floor(rgtau1/1000*fs));

[phi12_L, i_tau] = max(acfplot2_L(tg:end));
tau1_L = acfplot_L(i_tau+tg-1,1);
phi1_L = acfplot_L(i_tau+tg-1,2);

%W_phi(0)
for i = 1:tg;
    if acfplot_L(i,2) < 0.5
        tg2 = i;,break,
    end
end
pls = acfplot_L(tg2-1:tg2,:);
width_L = interp1(pls(:,2),pls(:,1),0.5,'pchip');

% --------------------------------------------------------------------
function [tau1_R phi1_R width_R] = acffactor1_R(tmax, fs, rgtau1, acfplot_R)

% calculation of tau_1 and phi_1
for i = 1:floor(rgtau1/1000*fs);
    if acfplot_R(i,2) < 0
        tg = i;,break,
    end
end

for i = 1:floor(rgtau1/1000*fs);
    tc(i,:) = exp(-acfplot_R(i,1)/0.01);
end
acfplot2_R = acfplot_R(1:floor(rgtau1/1000*fs),2).*tc(1:floor(rgtau1/1000*fs));

[phi12_R, i_tau] = max(acfplot2_R(tg:end));
tau1_R = acfplot_R(i_tau+tg-1,1);
phi1_R = acfplot_R(i_tau+tg-1,2);

%W_phi(0)
for i = 1:tg;
    if acfplot_R(i,2) < 0.5
        tg2 = i;,break,
    end
end
pls = acfplot_R(tg2-1:tg2,:);
width_R = interp1(pls(:,2),pls(:,1),0.5,'pchip');

% --------------------------------------------------------------------
function [taue_L corr_L index_taue_L reg0_L r_taue] = acffactor2_L(tmax, fs, acfplot_L, phi1_L, lgacfplot_L)

if phi1_L <= 0.1;
    acfplot2_L = flipud(acfplot_L(1:floor(fs*0.03),:));
    
    for i = 1:length(acfplot2_L);
        if acfplot2_L(i,2) > 0.1;
            taue_L = 1000*acfplot2_L(i,1);
            break
        end
    end
    
    corr_L = 0;
    index_taue_L = 0;
    reg0_L = 0;
    r_taue_L = 0;
    r_taue = zeros(10,2);
    
else
    
    % calculation of taue_L (every 3 ms)
    for i = 2:10
        r_taue_03 = lgacfplot_L(floor((i-1)*0.003*fs)+1:floor(i*0.003*fs),[1 2]);
        [p_taue_03(i,2), i_03] = max(r_taue_03(:,2));
        p_taue_03(i,1) = r_taue_03(i_03,1);
    end

    index_taue_03 = find(p_taue_03(:,2));
    r_03 = p_taue_03(index_taue_03,:);
    reg0_03 = polyfit(r_03(:,1),r_03(:,2),1);
    taue_03 = 1000*(-10-reg0_03(1,2))/reg0_03(1,1);
    cor_03 = -corrcoef(r_03(:,1),r_03(:,2));
    
    if size(cor_03,1) == 1,
        cor_03(1,2) = 0.0;
    else
        cor_03 = cor_03;
    end
    
    % calculation of taue_L (every 4 ms)
    for i = 2:10
        r_taue_04 = lgacfplot_L(floor((i-1)*0.004*fs)+1:floor(i*0.004*fs),[1 2]);
        [p_taue_04(i,2), i_04] = max(r_taue_04(:,2));
        p_taue_04(i,1) = r_taue_04(i_04,1);
    end

    index_taue_04 = find(p_taue_04(:,2));
    r_04 = p_taue_04(index_taue_04,:);
    reg0_04 = polyfit(r_04(:,1),r_04(:,2),1);
    taue_04 = 1000*(-10-reg0_04(1,2))/reg0_04(1,1);
    cor_04 = -corrcoef(r_04(:,1),r_04(:,2));
    
    if size(cor_04,1) == 1,
        cor_04(1,2) = 0.0;
    else
        cor_04 = cor_04;
    end
    
    % calculation of taue_L (every 5 ms)
    for i = 2:10
        r_taue_05 = lgacfplot_L(floor((i-1)*0.005*fs)+1:floor(i*0.005*fs),[1 2]);
        [p_taue_05(i,2), i_05] = max(r_taue_05(:,2));
        p_taue_05(i,1) = r_taue_05(i_05,1);
    end

    index_taue_05 = find(p_taue_05(:,2));
    r_05 = p_taue_05(index_taue_05,:);
    reg0_05 = polyfit(r_05(:,1),r_05(:,2),1);
    taue_05 = 1000*(-10-reg0_05(1,2))/reg0_05(1,1);
    cor_05 = -corrcoef(r_05(:,1),r_05(:,2));
    
    if size(cor_05,1) == 1,
        cor_05(1,2) = 0.0;
    else
        cor_05 = cor_05;
    end
    
    % calculation of taue_L (every 6 ms)
    for i = 2:10
        r_taue_06 = lgacfplot_L(floor((i-1)*0.006*fs)+1:floor(i*0.006*fs),[1 2]);
        [p_taue_06(i,2), i_06] = max(r_taue_06(:,2));
        p_taue_06(i,1) = r_taue_06(i_06,1);
    end

    index_taue_06 = find(p_taue_06(:,2));
    r_06 = p_taue_06(index_taue_06,:);
    reg0_06 = polyfit(r_06(:,1),r_06(:,2),1);
    taue_06 = 1000*(-10-reg0_06(1,2))/reg0_06(1,1);
    cor_06 = -corrcoef(r_06(:,1),r_06(:,2));
    
    if size(cor_06,1) == 1,
        cor_06(1,2) = 0.0;
    else
        cor_06 = cor_06;
    end
    
    % calculation of taue_L (every 7 ms)
    for i = 2:10
        r_taue_07 = lgacfplot_L(floor((i-1)*0.007*fs)+1:floor(i*0.007*fs),[1 2]);
        [p_taue_07(i,2), i_07] = max(r_taue_07(:,2));
        p_taue_07(i,1) = r_taue_07(i_07,1);
    end

    index_taue_07 = find(p_taue_07(:,2));
    r_07 = p_taue_07(index_taue_07,:);
    reg0_07 = polyfit(r_07(:,1),r_07(:,2),1);
    taue_07 = 1000*(-10-reg0_07(1,2))/reg0_07(1,1);
    cor_07 = -corrcoef(r_07(:,1),r_07(:,2));
    
    if size(cor_07,1) == 1,
        cor_07(1,2) = 0.0;
    else
        cor_07 = cor_07;
    end
    
    % calculation of taue_L (every 8 ms)
    for i = 2:10
        r_taue_08 = lgacfplot_L(floor((i-1)*0.008*fs)+1:floor(i*0.008*fs),[1 2]);
        [p_taue_08(i,2), i_08] = max(r_taue_08(:,2));
        p_taue_08(i,1) = r_taue_08(i_08,1);
    end

    index_taue_08 = find(p_taue_08(:,2));
    r_08 = p_taue_08(index_taue_08,:);
    reg0_08 = polyfit(r_08(:,1),r_08(:,2),1);
    taue_08 = 1000*(-10-reg0_08(1,2))/reg0_08(1,1);
    cor_08 = -corrcoef(r_08(:,1),r_08(:,2));
    
    if size(cor_08,1) == 1,
        cor_08(1,2) = 0.0;
    else
        cor_08 = cor_08;
    end
    
    % calculation of taue_L (every 9 ms)
    for i = 2:10
        r_taue_09 = lgacfplot_L(floor((i-1)*0.009*fs)+1:floor(i*0.009*fs),[1 2]);
        [p_taue_09(i,2), i_09] = max(r_taue_09(:,2));
        p_taue_09(i,1) = r_taue_09(i_09,1);
    end

    index_taue_09 = find(p_taue_09(:,2));
    r_09 = p_taue_09(index_taue_09,:);
    reg0_09 = polyfit(r_09(:,1),r_09(:,2),1);
    taue_09 = 1000*(-10-reg0_09(1,2))/reg0_09(1,1);
    cor_09 = -corrcoef(r_09(:,1),r_09(:,2));
    
    if size(cor_09,1) == 1,
        cor_09(1,2) = 0.0;
    else
        cor_09 = cor_09;
    end
    
    % calculation of taue_L (every 10 ms)
    for i = 2:10
        r_taue_10 = lgacfplot_L(floor((i-1)*0.010*fs)+1:floor(i*0.010*fs),[1 2]);
        [p_taue_10(i,2), i_10] = max(r_taue_10(:,2));
        p_taue_10(i,1) = r_taue_10(i_10,1);
    end

    index_taue_10 = find(p_taue_10(:,2));
    r_10 = p_taue_10(index_taue_10,:);
    reg0_10 = polyfit(r_10(:,1),r_10(:,2),1);
    taue_10 = 1000*(-10-reg0_10(1,2))/reg0_10(1,1);
    cor_10 = -corrcoef(r_10(:,1),r_10(:,2));
    
    if size(cor_10,1) == 1,
        cor_10(1,2) = 0.0;
    else
        cor_10 = cor_10;
    end
    
    set_taue =  [taue_03 taue_04 taue_05 taue_06 taue_07 taue_08 taue_09 taue_10
        cor_03(1,2) cor_04(1,2) cor_05(1,2) cor_06(1,2) cor_07(1,2) cor_08(1,2) cor_09(1,2) cor_10(1,2)];
    set_taue = set_taue';
    [corr_L, index_taue_L] = max(set_taue(:,2));
    taue_L = set_taue(index_taue_L,1);
    
    r_taue0 = zeros(10,16);
    r_taue0(1:length(r_03),[1 2]) = r_03;
    r_taue0(1:length(r_04),[3 4]) = r_04;
    r_taue0(1:length(r_05),[5 6]) = r_05;
    r_taue0(1:length(r_06),[7 8]) = r_06;
    r_taue0(1:length(r_07),[9 10]) = r_07;
    r_taue0(1:length(r_08),[11 12]) = r_08;
    r_taue0(1:length(r_09),[13 14]) = r_09;
    r_taue0(1:length(r_10),[15 16]) = r_10;
    r_taue = r_taue0(:,[index_taue_L*2-1 index_taue_L*2]);
    
    reg0b = [reg0_03; reg0_04; reg0_05; reg0_06; reg0_07; reg0_08; reg0_09; reg0_10];
    reg0_L = reg0b(index_taue_L,:);
    index_taue_L = index_taue_L+2;
    
    if corr_L < 0.85
        % calculation of taue (every 20 ms)
        for i = 1:8
            r_taue_20 = lgacfplot_L(floor((i-1)*0.020*fs)+1:floor(i*0.020*fs),[1 2]);
            [p_taue_20(i,2), i_20] = max(r_taue_20(:,2));
            p_taue_20(i,1) = r_taue_20(i_20,1);
        end
        for i = 2:4
            r_20(i,:) = p_taue_20(i,:);
        end
        
        if length(p_taue_20) > 4
            for i = 5:length(p_taue_20)
                if p_taue_20(i,2) < 10*log10(phi1_L)-5; break
                end
                r_20(i,:) = p_taue_20(i,:);
            end
        end
        
        index_taue_20 = find(r_20(:,2));
        r_20 = r_20(index_taue_20,:);
        num_20 = length(r_20);
        reg0_20 = polyfit(r_20(:,1),r_20(:,2),1);
        taue_20 = 1000*(-10-reg0_20(1,2))/reg0_20(1,1);
        cor_20 = -corrcoef(r_20(:,1),r_20(:,2));
        
        if size(cor_20,1) == 1,
            cor_20(1,2) = 0.0;
        else
            cor_20(1,2) = cor_20(1,2);
        end
        
        if cor_20(1,2) > corr_L,
            index_taue_L = 20;
            corr_L = cor_20(1,2);
            taue_L = taue_20;
            r_taue = zeros(10,2);
            r_taue(1:length(r_20),:) = r_20;
            reg0 = reg0_20;
        end
    end
    
    if corr_L < 0.85
        % calculation of taue (every 30 ms)
        for i = 1:5
            r_taue_30 = lgacfplot_L(floor((i-1)*0.030*fs)+1:floor(i*0.030*fs),[1 2]);
            [p_taue_30(i,2), i_30] = max(r_taue_30(:,2));
            p_taue_30(i,1) = r_taue_30(i_30,1);
        end
        for i = 2:4
            r_30(i,:) = p_taue_30(i,:);
        end
        
        if length(p_taue_30) > 4
            for i = 5:length(p_taue_30)
                if p_taue_30(i,2) < 10*log10(phi1_L)-5; break
                end
                r_30(i,:) = p_taue_30(i,:);
            end
        end
        num_30 = length(r_30);
        index_taue_30 = find(r_30(:,2));
        r_30 = r_30(index_taue_30,:);
        reg0_30 = polyfit(r_30(:,1),r_30(:,2),1);
        taue_30 = 1000*(-10-reg0_30(1,2))/reg0_30(1,1);
        cor_30 = -corrcoef(r_30(:,1),r_30(:,2));
        
        if size(cor_30,1) == 1,
            cor_30(1,2) = 0.0;
        else
            cor_30 = cor_30;
        end
        
        if corr_L < cor_30(1,2)
            corr_L = cor_30(1,2);
            taue_L = taue_30;
            r_taue = zeros(10,2);
            r_taue(1:length(r_30),:) = r_30;
            reg0 = reg0_30;
            index_taue_L = 30;
        end
    end
    
    if corr_L < 0.85
        % calculation of taue (every 40 ms)
        for i = 1:4
            r_taue_40 = lgacfplot_L(floor((i-1)*0.040*fs)+1:floor(i*0.040*fs),[1 2]);
            [p_taue_40(i,2), i_40] = max(r_taue_40(:,2));
            p_taue_40(i,1) = r_taue_40(i_40,1);
        end
        for i = 2:4
            r_40(i,:) = p_taue_40(i,:);
        end
        if length(p_taue_40) > 4
            for i = 5:length(p_taue_40)
                if p_taue_40(i,2) < 10*log10(phi1_L)-5; break
                end
                r_40(i,:) = p_taue_40(i,:);
            end
        end
        num_40 = length(r_40);
        index_taue_40 = find(r_40(:,2));
        r_40 = r_40(index_taue_40,:);
        reg0_40 = polyfit(r_40(:,1),r_40(:,2),1);
        taue_40 = 1000*(-10-reg0_40(1,2))/reg0_40(1,1);
        cor_40 = -corrcoef(r_40(:,1),r_40(:,2));
        
        if size(cor_40,1) == 1,
            cor_40(1,2) = 0.0;
        else
            cor_40 = cor_40;
        end
        
        if corr_L < cor_40(1,2)
            corr_L = cor_40(1,2);
            taue_L = taue_40;
            r_taue = zeros(10,2);
            r_taue(1:length(r_40),:) = r_40;
            reg0 = reg0_40;
            index_taue_L = 40;
        end
    end
    
        if corr_L < 0.85
            % calculation of taue (every 50 ms)
            for i = 1:4
                r_taue_50 = lgacfplot_L(floor((i-1)*0.050*fs)+1:floor(i*0.050*fs),[1 2]);
                [p_taue_50(i,2), i_50] = max(r_taue_50(:,2));
                p_taue_50(i,1) = r_taue_50(i_50,1);
            end
            for i = 2:4
                r_50(i,:) = p_taue_50(i,:);
            end
    
            num_50 = length(r_50);
            index_taue_50 = find(r_50(:,2));
            r_50 = r_50(index_taue_50,:);
            reg0_50 = polyfit(r_50(:,1),r_50(:,2),1);
            taue_50 = 1000*(-10-reg0_50(1,2))/reg0_50(1,1);
            cor_50 = -corrcoef(r_50(:,1),r_50(:,2));
    
            if size(cor_50,1) == 1,
                cor_50(1,2) = 0.0;
            else
                cor_50 = cor_50;
            end
    
            if corr_L < cor_50(1,2)
                corr_L = cor_50(1,2);
                taue_L = taue_50;
                r_taue = zeros(10,2);
                r_taue(1:length(r_50),:) = r_50;
                reg0 = reg0_50;
                index_taue_L = 50;
            end
        end
    
    if taue_L < 0;
        acfplot2_L = flipud(acfplot_L(1:floor(fs*0.03),:));
        
        for i = 1:length(acfplot2_L);
            if acfplot2_L(i,2) > 0.1;
                taue_L = 1000*acfplot2_L(i,1);
                break
            end
        end
        
        corr_L = 0;
        index_taue_L = 0;
        reg0_L = 0;
        r_taue_L = 0;
        
    end
    
end

% --------------------------------------------------------------------
function [taue_R corr_R index_taue_R reg0_R r_taue] = acffactor2_R(tmax, fs, acfplot_R, phi1_R, lgacfplot_R)

if phi1_R <= 0.1;
    acfplot2_R = flipud(acfplot_R(1:floor(fs*0.03),:));
    
    for i = 1:length(acfplot2_R);
        if acfplot2_R(i,2) > 0.1;
            taue_R = 1000*acfplot2_R(i,1);
            break
        end
    end
    
    corr_R = 0;
    index_taue_R = 0;
    reg0_R = 0;
    r_taue_R = 0;
    r_taue = zeros(10,2);
    
else
    
    % calculation of taue_R (every 3 ms)
    for i = 2:10
        r_taue_03 = lgacfplot_R(floor((i-1)*0.003*fs)+1:floor(i*0.003*fs),[1 2]);
        [p_taue_03(i,2), i_03] = max(r_taue_03(:,2));
        p_taue_03(i,1) = r_taue_03(i_03,1);
    end

    index_taue_03 = find(p_taue_03(:,2));
    r_03 = p_taue_03(index_taue_03,:);
    reg0_03 = polyfit(r_03(:,1),r_03(:,2),1);
    taue_03 = 1000*(-10-reg0_03(1,2))/reg0_03(1,1);
    cor_03 = -corrcoef(r_03(:,1),r_03(:,2));
    
    if size(cor_03,1) == 1,
        cor_03(1,2) = 0.0;
    else
        cor_03 = cor_03;
    end
    
    % calculation of taue_R (every 4 ms)
    for i = 2:10
        r_taue_04 = lgacfplot_R(floor((i-1)*0.004*fs)+1:floor(i*0.004*fs),[1 2]);
        [p_taue_04(i,2), i_04] = max(r_taue_04(:,2));
        p_taue_04(i,1) = r_taue_04(i_04,1);
    end

    index_taue_04 = find(p_taue_04(:,2));
    r_04 = p_taue_04(index_taue_04,:);
    reg0_04 = polyfit(r_04(:,1),r_04(:,2),1);
    taue_04 = 1000*(-10-reg0_04(1,2))/reg0_04(1,1);
    cor_04 = -corrcoef(r_04(:,1),r_04(:,2));
    
    if size(cor_04,1) == 1,
        cor_04(1,2) = 0.0;
    else
        cor_04 = cor_04;
    end
    
    % calculation of taue_R (every 5 ms)
    for i = 2:10
        r_taue_05 = lgacfplot_R(floor((i-1)*0.005*fs)+1:floor(i*0.005*fs),[1 2]);
        [p_taue_05(i,2), i_05] = max(r_taue_05(:,2));
        p_taue_05(i,1) = r_taue_05(i_05,1);
    end

    index_taue_05 = find(p_taue_05(:,2));
    r_05 = p_taue_05(index_taue_05,:);
    reg0_05 = polyfit(r_05(:,1),r_05(:,2),1);
    taue_05 = 1000*(-10-reg0_05(1,2))/reg0_05(1,1);
    cor_05 = -corrcoef(r_05(:,1),r_05(:,2));
    
    if size(cor_05,1) == 1,
        cor_05(1,2) = 0.0;
    else
        cor_05 = cor_05;
    end
    
    % calculation of taue_R (every 6 ms)
    for i = 2:10
        r_taue_06 = lgacfplot_R(floor((i-1)*0.006*fs)+1:floor(i*0.006*fs),[1 2]);
        [p_taue_06(i,2), i_06] = max(r_taue_06(:,2));
        p_taue_06(i,1) = r_taue_06(i_06,1);
    end

    index_taue_06 = find(p_taue_06(:,2));
    r_06 = p_taue_06(index_taue_06,:);
    reg0_06 = polyfit(r_06(:,1),r_06(:,2),1);
    taue_06 = 1000*(-10-reg0_06(1,2))/reg0_06(1,1);
    cor_06 = -corrcoef(r_06(:,1),r_06(:,2));
    
    if size(cor_06,1) == 1,
        cor_06(1,2) = 0.0;
    else
        cor_06 = cor_06;
    end
    
    % calculation of taue_R (every 7 ms)
    for i = 2:10
        r_taue_07 = lgacfplot_R(floor((i-1)*0.007*fs)+1:floor(i*0.007*fs),[1 2]);
        [p_taue_07(i,2), i_07] = max(r_taue_07(:,2));
        p_taue_07(i,1) = r_taue_07(i_07,1);
    end

    index_taue_07 = find(p_taue_07(:,2));
    r_07 = p_taue_07(index_taue_07,:);
    reg0_07 = polyfit(r_07(:,1),r_07(:,2),1);
    taue_07 = 1000*(-10-reg0_07(1,2))/reg0_07(1,1);
    cor_07 = -corrcoef(r_07(:,1),r_07(:,2));
    
    if size(cor_07,1) == 1,
        cor_07(1,2) = 0.0;
    else
        cor_07 = cor_07;
    end
    
    % calculation of taue_R (every 8 ms)
    for i = 2:10
        r_taue_08 = lgacfplot_R(floor((i-1)*0.008*fs)+1:floor(i*0.008*fs),[1 2]);
        [p_taue_08(i,2), i_08] = max(r_taue_08(:,2));
        p_taue_08(i,1) = r_taue_08(i_08,1);
    end

    index_taue_08 = find(p_taue_08(:,2));
    r_08 = p_taue_08(index_taue_08,:);
    reg0_08 = polyfit(r_08(:,1),r_08(:,2),1);
    taue_08 = 1000*(-10-reg0_08(1,2))/reg0_08(1,1);
    cor_08 = -corrcoef(r_08(:,1),r_08(:,2));
    
    if size(cor_08,1) == 1,
        cor_08(1,2) = 0.0;
    else
        cor_08 = cor_08;
    end
    
    % calculation of taue_R (every 9 ms)
    for i = 2:10
        r_taue_09 = lgacfplot_R(floor((i-1)*0.009*fs)+1:floor(i*0.009*fs),[1 2]);
        [p_taue_09(i,2), i_09] = max(r_taue_09(:,2));
        p_taue_09(i,1) = r_taue_09(i_09,1);
    end

    index_taue_09 = find(p_taue_09(:,2));
    r_09 = p_taue_09(index_taue_09,:);
    reg0_09 = polyfit(r_09(:,1),r_09(:,2),1);
    taue_09 = 1000*(-10-reg0_09(1,2))/reg0_09(1,1);
    cor_09 = -corrcoef(r_09(:,1),r_09(:,2));
    
    if size(cor_09,1) == 1,
        cor_09(1,2) = 0.0;
    else
        cor_09 = cor_09;
    end
    
    % calculation of taue_R (every 10 ms)
    for i = 2:10
        r_taue_10 = lgacfplot_R(floor((i-1)*0.010*fs)+1:floor(i*0.010*fs),[1 2]);
        [p_taue_10(i,2), i_10] = max(r_taue_10(:,2));
        p_taue_10(i,1) = r_taue_10(i_10,1);
    end

    index_taue_10 = find(p_taue_10(:,2));
    r_10 = p_taue_10(index_taue_10,:);
    reg0_10 = polyfit(r_10(:,1),r_10(:,2),1);
    taue_10 = 1000*(-10-reg0_10(1,2))/reg0_10(1,1);
    cor_10 = -corrcoef(r_10(:,1),r_10(:,2));
    
    if size(cor_10,1) == 1,
        cor_10(1,2) = 0.0;
    else
        cor_10 = cor_10;
    end
    
    set_taue =  [taue_03 taue_04 taue_05 taue_06 taue_07 taue_08 taue_09 taue_10
        cor_03(1,2) cor_04(1,2) cor_05(1,2) cor_06(1,2) cor_07(1,2) cor_08(1,2) cor_09(1,2) cor_10(1,2)];
    set_taue = set_taue';
    [corr_R, index_taue_R] = max(set_taue(:,2));
    taue_R = set_taue(index_taue_R,1);
    
    r_taue0 = zeros(10,18);
    r_taue0(1:length(r_03),[3 4]) = r_03;
    r_taue0(1:length(r_04),[5 6]) = r_04;
    r_taue0(1:length(r_05),[7 8]) = r_05;
    r_taue0(1:length(r_06),[9 10]) = r_06;
    r_taue0(1:length(r_07),[11 12]) = r_07;
    r_taue0(1:length(r_08),[13 14]) = r_08;
    r_taue0(1:length(r_09),[15 16]) = r_09;
    r_taue0(1:length(r_10),[17 18]) = r_10;
    r_taue = r_taue0(:,[index_taue_R*2-1 index_taue_R*2]);
    
    reg0b = [reg0_03; reg0_04; reg0_05; reg0_06; reg0_07; reg0_08; reg0_09; reg0_10];
    reg0_R = reg0b(index_taue_R,:);
    index_taue_R = index_taue_R+2;
    
    if corr_R < 0.85
        % calculation of taue (every 20 ms)
        for i = 1:10
            r_taue_20 = lgacfplot_R(floor((i-1)*0.020*fs)+1:floor(i*0.020*fs),[1 2]);
            [p_taue_20(i,2), i_20] = max(r_taue_20(:,2));
            p_taue_20(i,1) = r_taue_20(i_20,1);
        end
        for i = 2:4
            r_20(i,:) = p_taue_20(i,:);
        end
        
        if length(p_taue_20) > 4
            for i = 5:length(p_taue_20)
                if p_taue_20(i,2) < 10*log10(phi1_R)-5; break
                end
                r_20(i,:) = p_taue_20(i,:);
            end
        end
        
        index_taue_20 = find(r_20(:,2));
        r_20 = r_20(index_taue_20,:);
        num_20 = length(r_20);
        reg0_20 = polyfit(r_20(:,1),r_20(:,2),1);
        taue_20 = 1000*(-10-reg0_20(1,2))/reg0_20(1,1);
        cor_20 = -corrcoef(r_20(:,1),r_20(:,2));
        
        if size(cor_20,1) == 1,
            cor_20(1,2) = 0.0;
        else
            cor_20(1,2) = cor_20(1,2);
        end
        
        if cor_20(1,2) > corr_R,
            index_taue_R = 20;
            corr_R = cor_20(1,2);
            taue_R = taue_20;
            r_taue = zeros(10,2);
            r_taue(1:length(r_20),:) = r_20;
            reg0 = reg0_20;
        end
    end
    
    if corr_R < 0.85
        % calculation of taue (every 30 ms)
        for i = 1:6
            r_taue_30 = lgacfplot_R(floor((i-1)*0.030*fs)+1:floor(i*0.030*fs),[1 2]);
            [p_taue_30(i,2), i_30] = max(r_taue_30(:,2));
            p_taue_30(i,1) = r_taue_30(i_30,1);
        end
        for i = 2:4
            r_30(i,:) = p_taue_30(i,:);
        end
        
        if length(p_taue_30) > 4
            for i = 5:length(p_taue_30)
                if p_taue_30(i,2) < 10*log10(phi1_R)-5; break
                end
                r_30(i,:) = p_taue_30(i,:);
            end
        end
        num_30 = length(r_30);
        index_taue_30 = find(r_30(:,2));
        r_30 = r_30(index_taue_30,:);
        reg0_30 = polyfit(r_30(:,1),r_30(:,2),1);
        taue_30 = 1000*(-10-reg0_30(1,2))/reg0_30(1,1);
        cor_30 = -corrcoef(r_30(:,1),r_30(:,2));
        
        if size(cor_30,1) == 1,
            cor_30(1,2) = 0.0;
        else
            cor_30 = cor_30;
        end
        
        if corr_R < cor_30(1,2)
            corr_R = cor_30(1,2);
            taue_R = taue_30;
            r_taue = zeros(10,2);
            r_taue(1:length(r_30),:) = r_30;
            reg0 = reg0_30;
            index_taue_R = 30;
        end
    end
    
    if corr_R < 0.85
        % calculation of taue (every 40 ms)
        for i = 1:5
            r_taue_40 = lgacfplot_R(floor((i-1)*0.040*fs)+1:floor(i*0.040*fs),[1 2]);
            [p_taue_40(i,2), i_40] = max(r_taue_40(:,2));
            p_taue_40(i,1) = r_taue_40(i_40,1);
        end
        for i = 2:4
            r_40(i,:) = p_taue_40(i,:);
        end
        if length(p_taue_40) > 4
            for i = 5:length(p_taue_40)
                if p_taue_40(i,2) < 10*log10(phi1_R)-5; break
                end
                r_40(i,:) = p_taue_40(i,:);
            end
        end
        num_40 = length(r_40);
        index_taue_40 = find(r_40(:,2));
        r_40 = r_40(index_taue_40,:);
        reg0_40 = polyfit(r_40(:,1),r_40(:,2),1);
        taue_40 = 1000*(-10-reg0_40(1,2))/reg0_40(1,1);
        cor_40 = -corrcoef(r_40(:,1),r_40(:,2));
        
        if size(cor_40,1) == 1,
            cor_40(1,2) = 0.0;
        else
            cor_40 = cor_40;
        end
        
        if corr_R < cor_40(1,2)
            corr_R = cor_40(1,2);
            taue_R = taue_40;
            r_taue = zeros(10,2);
            r_taue(1:length(r_40),:) = r_40;
            reg0 = reg0_40;
            index_taue_R = 40;
        end
    end
    
    if corr_R < 0.85
        % calculation of taue (every 50 ms)
        for i = 1:4
            r_taue_50 = lgacfplot_R(floor((i-1)*0.050*fs)+1:floor(i*0.050*fs),[1 2]);
            [p_taue_50(i,2), i_50] = max(r_taue_50(:,2));
            p_taue_50(i,1) = r_taue_50(i_50,1);
        end
        for i = 2:4
            r_50(i,:) = p_taue_50(i,:);
        end
        
        num_50 = length(r_50);
        index_taue_50 = find(r_50(:,2));
        r_50 = r_50(index_taue_50,:);
        reg0_50 = polyfit(r_50(:,1),r_50(:,2),1);
        taue_50 = 1000*(-10-reg0_50(1,2))/reg0_50(1,1);
        cor_50 = -corrcoef(r_50(:,1),r_50(:,2));
        
        if size(cor_50,1) == 1,
            cor_50(1,2) = 0.0;
        else
            cor_50 = cor_50;
        end
        
        if corr_R < cor_50(1,2)
            corr_R = cor_50(1,2);
            taue_R = taue_50;
            r_taue = zeros(10,2);
            r_taue(1:length(r_50),:) = r_50;
            reg0 = reg0_50;
            index_taue_R = 50;
        end
    end
    
    if taue_R < 0;
        acfplot2_R = flipud(acfplot_R(1:floor(fs*0.03),:));
        
        for i = 1:length(acfplot2_R);
            if acfplot2_R(i,2) > 0.1;
                taue_R = 1000*acfplot2_R(i,1);
                break
            end
        end
        
        corr_R = 0;
        index_taue_R = 0;
        reg0_R = 0;
        r_taue_R = 0;
        
    end
    
end

% --------------------------------------------------------------------
function [IACC tau wiacc wpls wmns] = iacffactor(fs, ccfplot)

%IACC
tau_center = ceil(length(ccfplot)/2);
t1_lim = ceil(tau_center - 0.001*fs);
t2_lim = floor(tau_center + 0.001*fs);
range_iacf = ccfplot(t1_lim:t2_lim,:);
[IACC,index] = max(range_iacf(:,2));

%tau_IACC
tau  = ccfplot(t1_lim+index-1,1);

%W_IACC
dlt = 0.1;
mns0 = ccfplot(1:t1_lim+index-1,:);
mns = flipdim(mns0,1);
for i = 1:length(mns)
    mns2(i,:) = mns(i,:);
    if mns(i,2) < (1-dlt)*IACC
        mns2(i+1,:) = mns(i+1,:);
        break
    end
end

pls = ccfplot(t1_lim+index-1:end,:);
for i=1:length(pls)
    pls2(i,:) = pls(i,:);
    if pls(i,2) < (1-dlt)*IACC
        pls2(i+1,:) = pls(i+1,:);
        break
    end
end

wpls = interp1(mns2(1:length(mns2),2),mns2(1:length(mns2),1),(1-0.1)*IACC,'pchip');
wmns = interp1(pls2(1:length(pls2),2),pls2(1:length(pls2),1),(1-0.1)*IACC,'pchip');
wiacc = abs(wpls-wmns);

% --------------------------------------------------------------------
function y = afilter(x,fs);
%    AFILTER  Design of a A-weighting filter.
%    y = afilter(x,fs) designs a digital A-weighting filter for
%    sampling frequency fs.
%    Warning: fs should normally be higher than 20 kHz. For example,
%    fs = 48000 yields a class 1-compliant filter.
%
%    Requires the Signal Processing Toolbox.


% Definition of analog A-weighting filter according to IEC/CD 1672.

[B,A] = adsgn(fs);
y = filter(B,A,x);

% --------------------------------------------------------------------
function [B,A] = adsgn(fs)
% ADSGN  Design of a A-weighting filter.
%    [B,A] = ADSGN(fs) designs a digital A-weighting filter for
%    sampling frequency fs. Usage: Y = FILTER(B,A,X).
%    Warning: fs should normally be higher than 20 kHz. For example,
%    fs = 48000 yields a class 1-compliant filter.
%
%    Requires the Signal Processing Toolbox.
%
%    See also ASPEC, CDSGN, CSPEC.

% Author: Christophe Couvreur, Faculte Polytechnique de Mons (Belgium)
%         couvreur@thor.fpms.ac.be
% Last modification: Aug. 20, 1997, 10:00am.

% References:
%    [1] IEC/CD 1672: Electroacoustics-Sound Level Meters, Nov. 1996.

% Definition of analog A-weighting filter according to IEC/CD 1672.
f1 = 20.598997;
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;
A1000 = 1.9997;
pi = 3.14159265358979;
NUMs = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
DENs = conv(conv(DENs,[1 2*pi*f3]),[1 2*pi*f2]);

% Use the bilinear transformation to get the digital filter.
[B,A] = bilinear(NUMs,DENs,fs);
