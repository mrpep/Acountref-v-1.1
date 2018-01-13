function [TauMins,DeltaTau] = CalculateAnechoics(AN1,AN2,AN3,AN4,TauMinsAnechoics,fs)
    TauMins(1) = TauEMin(AN1,fs);
    TauMins(2) = TauEMin(AN2,fs);
    TauMins(3) = TauEMin(AN3,fs);
    TauMins(4) = TauEMin(AN4,fs);

    DeltaTau = TauMins - TauMinsAnechoics;

end