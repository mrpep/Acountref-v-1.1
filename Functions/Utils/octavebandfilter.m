%Filtrado por bandas de octava. Tambien devuelve las frecuencias centrales
%normalizadas

function [y,F0Normalized] = octavebandfilter(x,Fs)
    
    F0 = 1000;
    h = fdesign.octave(1,'Class 0', 'N,F0',6,F0,Fs);
    F0 = validfrequencies(h);
    F0Normalized = [32,63,125,250,500,1000,2000,4000,8000,16000];

    for i = 1:length(F0)
        y{i} = zeros(length(x),1);
        h.F0 = F0(i);
        H(i) = design(h,'butter');
        y{i} = filter(H(i),x);
        Fe(i) = Fs;        
    end
    
end