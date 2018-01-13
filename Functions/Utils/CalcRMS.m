%Devuelve el rms del archivo de calibracion x con frecuencia de muestreo fs
%usando tiempo de integracion int_time en segundos.

function [calMS] = CalcRMS(x,fs,int_time)

    calsq = x.^2;
    window = round(fs*int_time);
    k=1;
    for i=1:window:length(calsq)-window
        ms(k)=sum(calsq(i:i+window))/window;  
        k=k+1;
    end
    calMS = mean(ms);

end