function [LFs] = CalculateSF(IR,fs,octave,NC,threshold)

    IRBeginL = FindBeginIR(IR(:,1),threshold);
    IRBeginR = FindBeginIR(IR(:,2),threshold);
    if IRBeginL<IRBeginR
        IRBegin = IRBeginL;
    else
        IRBegin = IRBeginR;
    end
    IRL = IR(IRBegin:end,1);
    IRR = IR(IRBegin:end,2);
    
    [LFearly,LFlate] = LateralFraction(IRL,IRR,fs,octave);
    LFs = vertcat(LFearly,LFlate);
    
end

