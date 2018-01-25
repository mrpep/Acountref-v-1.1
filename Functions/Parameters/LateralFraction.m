function [LFearly,LFLate] = LateralFraction(IR8,IRomni,fs,octave)
    %Ventaneo:
    
    IR8early = IR8(floor(0.01*fs):floor(0.08*fs));
    IRomniearly = IRomni(1:floor(0.08*fs));
    IR8late = IR8(floor(0.08*fs):end);
    IRomnilate = IRomni(1:end);
    
    %Globales:
    LFearly(1) = 10*log10(trapz(IR8early.^2)/trapz(IRomniearly.^2));
    LFLate(1) = 10*log10(trapz(IR8late.^2)/trapz(IRomnilate.^2));
    
    %Subjetivos:
    
    IR810ms = integratewindow(IR8',fs,0.01);
    IRomni10ms = integratewindow(IRomni',fs,0.01);
    
    IR810early = IR810ms(floor(0.01*fs):floor(0.08*fs));
    IRomni10early = IRomni10ms(1:floor(0.08*fs));
    IR810late = IR810ms(floor(0.08*fs):end);
    IRomni10late = IRomni10ms(1:end);
    
    LFearly(2) = 10*log10(trapz(IR810early.^2)/trapz(IRomni10early.^2));
    LFLate(2) = 10*log10(trapz(IR810late.^2)/trapz(IRomni10late.^2));
    
    IR8100ms = integratewindow(IR8',fs,0.1);
    IRomni100ms = integratewindow(IRomni',fs,0.1);
    
    IR8100early = IR8100ms(floor(0.01*fs):floor(0.08*fs));
    IRomni100early = IRomni100ms(1:floor(0.08*fs));
    IR8100late = IR8100ms(floor(0.08*fs):end);
    IRomni100late = IRomni100ms(1:end);
    
    LFearly(3) = 10*log10(trapz(IR8100early.^2)/trapz(IRomni100early.^2));
    LFLate(3) = 10*log10(trapz(IR8100late.^2)/trapz(IRomni100late.^2));
    
    IR8350ms = integratewindow(IR8',fs,0.35);
    IRomni350ms = integratewindow(IRomni',fs,0.35);
    
    IR8350early = IR8350ms(floor(0.01*fs):floor(0.08*fs));
    IRomni350early = IRomni350ms(1:floor(0.08*fs));
    IR8350late = IR8350ms(floor(0.08*fs):end);
    IRomni350late = IRomni350ms(1:end);
    
    LFearly(4) = 10*log10(trapz(IR8350early.^2)/trapz(IRomni350early.^2));
    LFLate(4) = 10*log10(trapz(IR8350late.^2)/trapz(IRomni350late.^2));
    
    
    %Filtro las ventanas para evitar group delay:
    [IR8earlyfiltered,Fs,fc] = FilterThirdOctave(IR8early,fs,octave);
    [IRomniearlyfiltered,~,~] = FilterThirdOctave(IRomni,fs,octave);
    [IR8latefiltered,~,~] = FilterThirdOctave(IR8late,fs,octave);
    [IRomnilatefiltered,~,~] = FilterThirdOctave(IRomnilate,fs,octave);
    
    %Computo por bandas:
    for i=1:length(IR8earlyfiltered)
        LFearly(i+4) = 10*log10(trapz(IR8earlyfiltered{i}.^2)/trapz(IRomniearlyfiltered{i}.^2));
        LFLate(i+4) = 10*log10(trapz(IR8latefiltered{i}.^2)/trapz(IRomnilatefiltered{i}.^2));
    end
    
end