%Esta función computa la división entre dos integrales de la presion al
%cuadrado. a,b,c y d son los limites de integración ([a,b] para integral
%del numerador, [c,d] para integral del denominador). -1 es infinito
%(ultima muestra de la señal)

function [energyratio] = energyratioparam(IR,IR10ms,IR100ms,IR350ms,fs,a,b,c,d,octave)

    asamples = round(fs*a/1000); bsamples = round(fs*b/1000); csamples = round(fs*c/1000); dsamples = round(fs*d/1000);
    bsamples10 = bsamples; bsamples100 = bsamples; bsamples350 = bsamples;
    dsamples10 = dsamples; dsamples100 = dsamples; dsamples350 = dsamples;

    if (a == 0)
        asamples = 1;
    end
    if (b == -1)
        
        bsamples = length(IR);
        bsamples10 = length(IR10ms);
        bsamples100 = length(IR100ms);
        bsamples350 = length(IR350ms);
        
    end
    if (c == 0)
        csamples = 1;
    end
    if (d == -1)
        
        dsamples = length(IR);
        dsamples10 = length(IR10ms);
        dsamples100 = length(IR100ms);
        dsamples350 = length(IR350ms);
        
    end
    if (dsamples > length(IR))
        dsamples = length(IR);
    end
    if (dsamples10 > length(IR10ms))
        dsamples10 = length(IR10ms);
    end
    if (dsamples100 > length(IR100ms))
        dsamples100 = length(IR100ms);
    end
    if (dsamples350 > length(IR350ms))
        dsamples350 = length(IR350ms);
    end
    
    if (bsamples > length(IR))
        bsamples = length(IR);
    end
    if (bsamples10 > length(IR10ms))
        bsamples10 = length(IR10ms);
    end
    if (bsamples100 > length(IR100ms))
        bsamples100 = length(IR100ms);
    end
    if (bsamples350 > length(IR350ms))
        bsamples350 = length(IR350ms);
    end
    
    %Primero ventaneo la señal
    num = IR(asamples:bsamples);
    den = IR(csamples:dsamples);
    intnum = trapz(num.^2);
    intden = trapz(den.^2);

    energyratio(1) = 10*log10(intnum/intden); %Valor global

    %Cálculo subjetivo:
    %10 ms:
    IR10num = IR10ms(asamples:bsamples10);
    IR10den = IR10ms(csamples:dsamples10);
    intnum10 = trapz(IR10num.^2);
    intden10 = trapz(IR10den.^2);
    energyratio(2) = 10*log10(intnum10/intden10);
    
    %100ms:
        
    IR100num = IR100ms(asamples:bsamples100);
    IR100den = IR100ms(csamples:dsamples100);
    intnum100 = trapz(IR100num.^2);
    intden100 = trapz(IR100den.^2);
    energyratio(3) = 10*log10(intnum100/intden100);
    
    %350 ms:
    
    IR350num = IR350ms(asamples:bsamples350);
    IR350den = IR350ms(csamples:dsamples350);
    intnum350 = trapz(IR350num.^2);
    intden350 = trapz(IR350den.^2);
    energyratio(4) = 10*log10(intnum350/intden350);
    
    %Filtro las ventanas para evitar group delay:
    [numFiltered,Fs,fc] = FilterThirdOctave(num,fs,octave);
    [denFiltered,~,~] = FilterThirdOctave(den,fs,octave);

    for i=1:length(numFiltered)
        intnum = trapz(numFiltered{i}.^2);
        intden = trapz(denFiltered{i}.^2);
        energyratio(4+i) = 10*log10(eps+intnum/intden); %Valor global
    end

end