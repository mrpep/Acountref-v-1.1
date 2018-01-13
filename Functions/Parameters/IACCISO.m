%Según ISO:

function [IACFearly,IACFlate,IACFall,IACCearly,IACClate,IACCall] = IACCISO(IRL,IRR,fs)
idx80 = fs*0.08;
idx500 = fs*0.5;
IACFearly = xcorr(IRL(1:idx80),IRR(1:idx80))/(trapz(IRL(1:idx80).^2)*trapz(IRR(1:idx80).^2)).^(0.5);
IACFlate = xcorr(IRL(idx80:idx500),IRR(idx80:idx500))/(trapz(IRL(idx80:idx500).^2)*trapz(IRR(idx80:idx500).^2)).^(0.5);
IACFall = xcorr(IRL(1:idx500),IRR(1:idx500))/(trapz(IRL(1:idx500).^2)*trapz(IRR(1:idx500).^2)).^(0.5);
IACCearly = max(IACFearly);
IACClate = max(IACFlate);
IACCall = max(IACFall);
end