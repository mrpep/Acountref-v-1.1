function [crosspoint,fmodel,ecm] = PepinoMethod4(t,y,fs,resolution)
    
    actualresolution = 0.2;
    %samplesiter = floor(actualresolution*fs);
    leftpoint = actualresolution;
    rightpoint = t(end)-actualresolution;
    steperror = 1;
    
while (actualresolution>resolution)
    
    crosspointgrid = leftpoint:actualresolution:rightpoint;
    ecm = zeros(1,length(crosspointgrid));
    
    for k=1:length(crosspointgrid)
        crosspointsample = find(t<=crosspointgrid(k),1,'last');
        modeldecaycoef = polyfit(t(1:crosspointsample),y(1:crosspointsample),1);
        modeldecay = polyval(modeldecaycoef,t(1:crosspointsample));
        modelnoise = ones(1,length(y)-crosspointsample)*modeldecay(end);
        model = horzcat(modeldecay,modelnoise);
        ecm(k) = sum((model(1:steperror:end)-y(1:end)).^2)/length(y);
    end
    
    crosspointk = find(ecm==min(ecm));
    
    actualresolution = actualresolution/10;
    %samplesiter = floor(actualresolution*fs);
    leftpoint = crosspointgrid(crosspointk)-2*actualresolution;
    rightpoint = crosspointgrid(crosspointk) + 2*actualresolution;
    if (leftpoint<actualresolution)
        leftpoint = actualresolution;
    end
    if (rightpoint>length(y)-actualresolution)
        rightpoint = length(y)-actualresolution;
    end
end
 
    crosspoint = crosspointgrid(crosspointk);  
    crosspointsample = find(t<=crosspointgrid(crosspointk),1,'last');
    fmodeldecaycoef = polyfit(t(1:crosspointsample),y(1:crosspointsample),1);
    fmodeldecay = polyval(fmodeldecaycoef,t(1:crosspointsample));
    fmodelnoise = ones(1,length(y)-crosspointsample)*fmodeldecay(end);
    fmodel = horzcat(fmodeldecay,fmodelnoise);
   
    
end