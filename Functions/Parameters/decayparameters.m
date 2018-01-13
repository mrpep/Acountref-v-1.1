function [EDT,T20,T30] = decayparameters(schroeder,fs)

%Index for TR calculus
    i5 = find(schroeder > schroeder(1)-5,1,'last');
    i15 = find(schroeder > schroeder(1)-10,1,'last');
    i25 = find(schroeder > schroeder(1)-25,1,'last');
    i35 = find(schroeder > schroeder(1)-35,1,'last');
    
t = (0:length(schroeder)-1)/fs;
tedt = t(1:i15);
t20 = t(i5:i25);
t30 = t(i5:i35);

edtc = polyfit(tedt,schroeder(1:i15),1);
t20c = polyfit(t20,schroeder(i5:i25),1);
t30c = polyfit(t30,schroeder(i5:i35),1);

EDT = -60/edtc(1);
T20 = -60/t20c(1);
T30 = -60/t30c(1);

end