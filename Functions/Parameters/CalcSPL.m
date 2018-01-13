%Calcula a partir del RMS de la calibracion, el audio a analizar x con su
%frecuencia de muestreo fs, el SPL de referencia del audio de calibracion y
%el tiempo de integracion int_time, el vector SPLTime con el logger de Leqs
%ventana a ventana, y el valor Leq con el Leq de todo el audio.

function [Leq] = CalcSPL(x,fs,calMS,refSPL,int_time)
  
    xsq = x.^2;
    window = round(fs*int_time);
    k=1;
    for i=1:window:length(xsq)-window
        ms(k)=sum(xsq(i:i+window))/window;  
        SPLTime(k) = refSPL + 10*log10(ms(k)/calMS); 
        k=k+1;
    end
    Leq = 10*log10(sum(10.^(SPLTime/10))/length(SPLTime));
    
end