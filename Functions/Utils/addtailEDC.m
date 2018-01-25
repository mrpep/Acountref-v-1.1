function [CorrectedSchroeder] = addtailEDC(schroeder,fs)
tlate = (floor(0.7*length(schroeder)):floor(0.9*length(schroeder))-1)/fs;
coeflate = polyfit(tlate,schroeder(floor(0.7*length(schroeder))+1:floor(0.9*length(schroeder))),1);
toaddsamples = -ceil((35/coeflate(1))*fs);
ttoadd = (floor(0.9*length(schroeder)):length(schroeder)+toaddsamples)/fs;
compensation = polyval(coeflate,ttoadd);
CorrectedSchroeder = horzcat(schroeder(1:floor(0.9*length(schroeder))),compensation);

end