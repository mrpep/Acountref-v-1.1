function [crosspoint,fmodel,ecm] = PepinoMethodSlow(t,y,fs,resolution)
    
    y = 10*log10(y.^2 + eps());

    crosspointgrid = (resolution:resolution:t(end));
    ecm = zeros(1,length(crosspointgrid));

    for k=1:length(crosspointgrid)
        crosspointsample = find(t<=crosspointgrid(k),1,'last');
        modeldecaycoef = polyfit(t(1:crosspointsample),y(1:crosspointsample),1);
        modeldecay = polyval(modeldecaycoef,t(1:crosspointsample));
        modelnoise = ones(1,length(y)-crosspointsample)*modeldecay(end);
        model = horzcat(modeldecay,modelnoise);
        ecm(k) = sum((model-y).^2)/length(y);
    end

    crosspointk = find(ecm==min(ecm));   
    crosspoint = crosspointgrid(crosspointk);  
    crosspointsample = find(t<=crosspointgrid(crosspointk),1,'last');
    fmodeldecaycoef = polyfit(t(1:crosspointsample),y(1:crosspointsample),1);
    fmodeldecay = polyval(fmodeldecaycoef,t(1:crosspointsample));
    fmodelnoise = ones(1,length(y)-crosspointsample)*fmodeldecay(end);
    fmodel = horzcat(fmodeldecay,fmodelnoise);
    
end