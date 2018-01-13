%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intlinear.m
%
%[A,B] = intlinear(x,y)
%
% Does de linear regression (least mean squares) for given (x,y).
% Returns A and B, with meaning y = A + B*x.

function    [A,B] = intlinear(x,y)

mx = nanmean(x);        my = nanmean(y);
mx2 = nanmean(x.^2);    my2 = nanmean(y.^2);
mxy = nanmean(x.*y);
A = (mx2*my-mx*mxy)/(mx2-mx^2);
B = (mxy - (mx*my))/(mx2-mx^2);

end