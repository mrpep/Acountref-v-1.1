function [Parameters,fc,Fs] = CalculateParametersStereo(IR,AN1,AN2,AN3,AN4,fs,CalFactor,octave,NC,threshold)

    IRBeginL = FindBeginIR(IR(:,1),threshold);
    IRBeginR = FindBeginIR(IR(:,2),threshold);
    if IRBeginL<IRBeginR
        IRBegin = IRBeginL;
    else
        IRBegin = IRBeginR;
    end
    IRL = IR(IRBegin:end,1);
    IRR = IR(IRBegin:end,2);
    
    [IRLFiltered,~,~] = FilterThirdOctave(IRL,fs,octave); %Filtro que no es de Matlab
    [IRRFiltered,Fs,fc] = FilterThirdOctave(IRR,fs,octave); %Filtro que no es de Matlab
    
    %Globales:
    [IACFearly,IACFlate,IACFall,IACCearly,IACClate,IACCall] = IACCISO(IRL,IRR,fs);
    IACFearlys{1} = IACFearly;
    IACFlates{1} = IACFlate;
    IACFalls{1} = IACFall;
    IACCearlys(1) = IACCearly;
    IACClates(1) = IACClate;
    IACCalls(1) = IACCall;
    IACCAndo = CalculateAnechoicsStereo(AN1,AN2,AN3,AN4,fs);
    
    for i=1:length(IRLFiltered)
        
        [IACFearly,IACFlate,IACFall,IACCearly,IACClate,IACCall] = IACCISO(IRLFiltered{i},IRRFiltered{i},Fs(i));
        IACFearlys{i+1} = IACFearly;
        IACFlates{i+1} = IACFlate;
        IACFalls{i+1} = IACFall;
        IACCearlys(i+1) = IACCearly;
        IACClates(i+1) = IACClate;
        IACCalls(i+1) = IACCall;       
    end
        
    IACCAndos = zeros(4,length(IACCearlys));
    IACCAndos(1,1) = IACCAndo(1);
    IACCAndos(2,1) = IACCAndo(2);
    IACCAndos(3,1) = IACCAndo(3);
    IACCAndos(4,1) = IACCAndo(4);
    
    Parameters = vertcat(IACCearlys,IACClates,IACCalls,IACCAndos);
end