function [IACC] = IACCAndo(x,fs)

    RunningFactor = ACF(x,fs,1,0.1,0.2,10);
    IACC = mean(RunningFactor(:,6)); %Es realmente esto lo que se busca?
end

