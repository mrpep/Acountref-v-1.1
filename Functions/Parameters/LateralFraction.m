function [LFearly,LFLate] = LateralFraction(IR8,IRomni,fs,octave)
    %Ventaneo:
    
    IR8early = IR8(floor(0.01*fs):floor(0.08*fs));
    IRomniearly = IRomni(1:floor(0.08*fs));
    IR8late = IR8(floor(0.08*fs):end);
    IRomnilate = IRomni(1:end);
    
    %Globales:
    LFearly(1) = 10*log10(trapz(IR8early.^2)/trapz(IRomniearly.^2));
    LFLate(1) = 10*log10(trapz(IR8late.^2)/trapz(IRomnilate.^2));
    
    %Filtro las ventanas para evitar group delay:
    [IR8earlyfiltered,Fs,fc] = FilterThirdOctave(IR8early,fs,octave);
    [IRomniearlyfiltered,~,~] = FilterThirdOctave(IRomni,fs,octave);
    [IR8latefiltered,~,~] = FilterThirdOctave(IR8late,fs,octave);
    [IRomnilatefiltered,~,~] = FilterThirdOctave(IRomnilate,fs,octave);
    
    %Computo por bandas:
    for i=1:length(IR8earlyfiltered)
        LFearly(i+1) = 10*log10(trapz(IR8earlyfiltered{i}.^2)/trapz(IRomniearlyfiltered{i}.^2));
        LFLate(i+1) = 10*log10(trapz(IR8latefiltered{i}.^2)/trapz(IRomnilatefiltered{i}.^2));
    end
    
end