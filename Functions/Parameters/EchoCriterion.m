function [EchoCurve,maxecho,tmax] = EchoCriterion(x,fs,n,te)
tesamples = round(te*fs);


numsum = 0;
densum = 0;
for tau = 1:length(x)
    numsum = numsum + (tau/fs)*(abs(x(tau)))^n;
    densum = densum + (abs(x(tau)))^n;
    ts(tau) = numsum/densum;
end

for i=(tesamples+1):length(ts)
    EchoCurve(i) = (ts(i)-ts(i-tesamples))/(tesamples/fs);
end
maxecho = max(EchoCurve);
tmax = find(EchoCurve==max(EchoCurve),1);
tmax = tmax/fs;
