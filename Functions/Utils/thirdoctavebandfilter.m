%Filtrado por tercios de octava. Tambien devuelve las frecuencias centrales
%normalizadas

function [y,F0Normalized] = thirdoctavebandfilter(x,Fs)
    
    F0 = 1000;
    h = fdesign.octave(3,'Class 0', 'N,F0',8,F0,Fs);
    F0 = validfrequencies(h);
    F0Normalized = [25,32,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1260,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000];

    for i = 1:length(F0)
        y{i} = zeros(length(x),1);
        h.F0 = F0(i);
        H(i) = design(h,'butter');
        y{i} = filter(H(i),x);
        Fe(i) = Fs;        
    end
    
end