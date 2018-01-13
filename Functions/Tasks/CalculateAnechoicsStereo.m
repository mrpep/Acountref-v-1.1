function [IACC] = CalculateAnechoicsStereo(AN1,AN2,AN3,AN4,fs)
    IACC(1) = IACCAndo(AN1,fs);
    IACC(2) = IACCAndo(AN2,fs);
    IACC(3) = IACCAndo(AN3,fs);
    IACC(4) = IACCAndo(AN4,fs);
end

