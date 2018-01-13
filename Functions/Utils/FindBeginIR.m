function[begink] = FindBeginIR(x,threshold)
    x = 10*log10(x.^2);
    x = x-max(x);
    begink = find(x>=max(x)+threshold,1);
    %cutIR = x(begink:end);
    
end