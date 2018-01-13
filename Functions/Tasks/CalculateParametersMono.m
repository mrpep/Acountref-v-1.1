function [Parameters,EchoSpeechs,EchoMusics,fc,Fs] = CalculateParametersMono(IR,WN,fs,octave,NC,threshold,calRMS)
    %1) Encontrar comienzo IR, cortar
    
    IRBegin = FindBeginIR(IR,threshold);
    IRcut = IR(IRBegin:end);
        
    %2) Filtrar:
    [IRFiltered,Fs,fc] = FilterThirdOctave(IRcut,fs,octave); %Filtro que no es de Matlab
    [WNFiltered,FsW,fcW] = FilterThirdOctave(WN,fs,octave);
    t=(0:length(IRcut)-1)/fs;
    %3) Correccion ruido        
    if NC == 1
        crosspoint = lundeby(IRcut,fs);
    elseif NC == 2
        crosspoint = fs*PepinoMethodFast(t,IRcut,fs,0.01);
    elseif NC == 0
        crosspoint = length(IRcut);
    end
    % Calculo de globales (sin filtrar)
    
    IRglobal = IRcut(1:crosspoint);
    tic
    IR10ms = integratewindow(IRglobal,fs,0.01);
    ms10 = toc
    tic
    IR100ms = integratewindow(IRglobal,fs,0.1);
    ms100 = toc
    tic
    IR350ms = integratewindow(IRglobal,fs,0.35);
    ms350 = toc
    
    C50 = energyratioparam(IRglobal,IR10ms,IR100ms,IR350ms,fs,0,50,50,-1,octave);
    C80 = energyratioparam(IRglobal,IR10ms,IR100ms,IR350ms,fs,0,80,80,-1,octave);
    D50 = energyratioparam(IRglobal,IR10ms,IR100ms,IR350ms,fs,0,50,0,-1,octave);
    DR = energyratioparam(IRglobal,IR10ms,IR100ms,IR350ms,fs,0,3,3,-1,octave);
    StageSupportEarly1 = energyratioparam(IRglobal,IR10ms,IR100ms,IR350ms,fs,20,100,0,10,octave);
    StageSupportEarly2 = energyratioparam(IRglobal,IR10ms,IR100ms,IR350ms,fs,20,200,0,10,octave);
    StageSupportLate = energyratioparam(IRglobal,IR10ms,IR100ms,IR350ms,fs,100,1000,0,10,octave); 
    RunningReverberance = energyratioparam(IRglobal,IR10ms,IR100ms,IR350ms,fs,160,320,0,160,octave); 
    IRsq = IRglobal.^2;
    edc = 10*log10(fliplr(cumtrapz(IRsq(end:-1:1))));
    [EDT,T20,T30] = decayparameters(edc,fs);
    [EchoSpeech,maxechospeech,tmaxspeech] = EchoCriterion(IRglobal,fs,2/3,0.009);
    [EchoMusic,maxechomusic,tmaxmusic] = EchoCriterion(IRglobal,fs,1,0.014);
    EDTs(1) = EDT;
    T20s(1) = T20;
    T30s(1) = T30;
    EchoMusics{1} = EchoMusic;
    EchoMusicMaxs(1) = maxechomusic;
    EchoMusicTMaxs(1) = tmaxmusic;
    EchoSpeechs{1} = EchoSpeech;
    EchoSpeechMaxs(1) = maxechospeech;
    EchoSpeechTMaxs(1) = tmaxspeech;
    SPLs(1) = CalcSPL(WN,fs,calRMS,94,1);
    %Computo de STI y ALCONS%:
    [STIval, ALcons,~,~] = STI(IRcut,fs);
    STIs(1) = STIval;
    ALconss(1) = ALcons;

    %Cálculo subjetivo:
    %10 ms:
        
    IRsq10 = IR10ms.^2;
    edc10 = 10*log10(fliplr(cumtrapz(IRsq10(end:-1:1))));
    [EDT10,T2010,T3010] = decayparameters(edc10,fs);
    [EchoSpeech10,maxechospeech10,tmaxspeech10] = EchoCriterion(IR10ms,fs,2/3,0.009);
    [EchoMusic10,maxechomusic10,tmaxmusic10] = EchoCriterion(IR10ms,fs,1,0.014);
    EDTs(2) = EDT10;
    T20s(2) = T2010;
    T30s(2) = T3010;
    EchoMusics{2} = EchoMusic10;
    EchoMusicMaxs(2) = maxechomusic10;
    EchoMusicTMaxs(2) = tmaxmusic10;
    EchoSpeechs{2} = EchoSpeech10;
    EchoSpeechMaxs(2) = maxechospeech10;
    EchoSpeechTMaxs(2) = tmaxspeech10;
    SPLs(2) = 0;
    %Computo de STI y ALCONS%:
    [STIval10, ALcons10,~,~] = STI(IR10ms,fs);
    STIs(2) = STIval10;
    ALconss(2) = ALcons10;
    
    %100 ms:
        
    IRsq100 = IR100ms.^2;
    edc100 = 10*log10(fliplr(cumtrapz(IRsq100(end:-1:1))));
    [EDT100,T20100,T30100] = decayparameters(edc100,fs);
    [EchoSpeech100,maxechospeech100,tmaxspeech100] = EchoCriterion(IR100ms,fs,2/3,0.009);
    [EchoMusic100,maxechomusic100,tmaxmusic100] = EchoCriterion(IR100ms,fs,1,0.014);
    EDTs(3) = EDT100;
    T20s(3) = T20100;
    T30s(3) = T30100;
    EchoMusics{3} = EchoMusic100;
    EchoMusicMaxs(3) = maxechomusic100;
    EchoMusicTMaxs(3) = tmaxmusic100;
    EchoSpeechs{3} = EchoSpeech100;
    EchoSpeechMaxs(3) = maxechospeech100;
    EchoSpeechTMaxs(3) = tmaxspeech100;
    SPLs(3) = 0;
    %Computo de STI y ALCONS%:
    [STIval100, ALcons100,~,~] = STI(IR100ms',fs);
    STIs(3) = STIval100;
    ALconss(3) = ALcons100;
    
    %350 ms:
    
    
    IRsq350 = IR350ms.^2;
    edc350 = 10*log10(fliplr(cumtrapz(IRsq350(end:-1:1))));
    [EDT350,T20350,T30350] = decayparameters(edc350,fs);
    [EchoSpeech350,maxechospeech350,tmaxspeech350] = EchoCriterion(IR350ms,fs,2/3,0.009);
    [EchoMusic350,maxechomusic350,tmaxmusic350] = EchoCriterion(IR350ms,fs,1,0.014);
    EDTs(4) = EDT350;
    T20s(4) = T20350;
    T30s(4) = T30350;
    EchoMusics{4} = EchoMusic350;
    EchoMusicMaxs(4) = maxechomusic350;
    EchoMusicTMaxs(4) = tmaxmusic350;
    EchoSpeechs{4} = EchoSpeech350;
    EchoSpeechMaxs(4) = maxechospeech350;
    EchoSpeechTMaxs(4) = tmaxspeech350;
    SPLs(4) = 0;
    %Computo de STI y ALCONS%:
    [STIval350, ALcons350,~,~] = STI(IR350ms,fs);
    STIs(4) = STIval350;
    ALconss(4) = ALcons350;
            
    %Calculo por octavas
    for i=1:length(IRFiltered)
        
        IRBegin = FindBeginIR(IRFiltered{i},threshold);
        IRFiltered{i} = IRFiltered{i}(IRBegin:end);
        IRFilteredETC{i} = 10*log10(IRFiltered{i}.^2+eps());
        t=(0:length(IRFilteredETC{i})-1)/Fs(i); %Chequear si se puede sacar del for.
        
        %3) Correccion ruido        
        if NC == 1
            crosspoint = Fs(i)*lundeby(IRFilteredETC{i},fs);
        elseif NC == 2
            crosspoint = Fs(i)*PepinoMethodFast(t,IRFilteredETC{i},fs,0.01);
        elseif NC == 0
            crosspoint = length(IRFiltered{i});
        end
        IRFiltered{i} = IRFiltered{i}(1:crosspoint);
        
        %4) Cálculo parámetros

            SPLs(i+4) = CalcSPL(WNFiltered{i},FsW(i),calRMS,94,1);
            
            %b) Decay Rates:
            IRFilteredsq = IRFiltered{i}.^2;
            edc = 10*log10(fliplr(cumtrapz(IRFilteredsq(end:-1:1))));
            [EDT,T20,T30] = decayparameters(edc,Fs(i));
            EDTs(i+4) = EDT;
            T20s(i+4) = T20;
            T30s(i+4) = T30;
            
            %c) Intelligibility:
            STIs(i+4) = 0;
            ALconss(i+4) = 0;
            
            %d) Echo Criterias:
            [EchoSpeech,maxechospeech,tmaxspeech] = EchoCriterion(IRFiltered{i},Fs(i),2/3,0.009);
            [EchoMusic,maxechomusic,tmaxmusic] = EchoCriterion(IRFiltered{i},Fs(i),1,0.014);
            
            EchoMusics{i+4} = EchoMusic;
            EchoMusicMaxs(i+4) = maxechomusic;
            EchoMusicTMaxs(i+4) = tmaxmusic;
            EchoSpeechs{i+4} = EchoSpeech;
            EchoSpeechMaxs(i+4) = maxechospeech;
            EchoSpeechTMaxs(i+4) = tmaxspeech;
    end  
    
    Parameters = vertcat(C50,C80,D50,EDTs,T20s,T30s,STIs,ALconss,EchoMusicMaxs,EchoMusicTMaxs,EchoSpeechMaxs,EchoSpeechTMaxs,DR,StageSupportEarly1,StageSupportEarly2,StageSupportLate,RunningReverberance,SPLs);
         
end