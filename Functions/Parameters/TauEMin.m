function [TauMin] = TauEMin(x,fs)

    RunningFactor = ACF(x,fs,1,0.1,0.2,10);
    TauMin = prctile(RunningFactor(:,6),5); %percentil 95%
    
end