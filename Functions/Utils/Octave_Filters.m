function Filter3 = Octave_Filters(Fs)

%%%%%%%%%%%%%%%
%
% Filter3 = Octave_Filters(Fs)
% 
% FUNCTION
%     This function designs IIR filters for each third octave bands with
%     the following center frequencies
%     Fc = [31.5 63 125 250 500 1000 2000 4000 8000 16000];
%
%   INPUT
%   Fs : sampling frequency in Hz
%        (NOTE: low frequencies filters may have a lower sampling frequency)
%
%
%   OUTPUT
%   Filter3: array of 10 structures which contains the filters
%            Filter3(i).FS: sampling freq. of filter for band n. i
%            Filter3(i).A : denomminator coefficients of IIR filter for band n. i
%            Filter3(i).B : numerator coefficients of IIR filter for band n. i
%            Filter3(i).Fc: center frequency of band n. i
%
% NOTE
%       this function requires OCTDSGN function from 3rd oct. toolbox of
%       Christophe Couvreur (downloadable from MatlabCentral)
%
%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - 2009
%%%%%%%%%%%%%%%%%%%

%%	Begin function
Fc = [31.5 63 125 250 500 1000 2000 4000 8000 ];

%%	Filters with fc < 220 will be resampled 
FiltOrd=4;
Filter3 = [];

%%	Filter resampled for fc = 31.5 Hz.
ink=1;
q=8; 	
FsNew=100*floor(Fs/q/100);
[B,A] = octdsgn(Fc(ink),FsNew,FiltOrd);
   Filter3(ink).FS = FsNew;
   Filter3(ink).A = A;
   Filter3(ink).B = B;
   Filter3(ink).Fc = Fc(ink);  

%%	Filter resampled for fc = 63 Hz.
ink=2;
q=4; 	
FsNew=100*floor(Fs/q/100);
[B,A] = octdsgn(Fc(ink),FsNew,FiltOrd);
   Filter3(ink).FS = FsNew;
   Filter3(ink).A = A;
   Filter3(ink).B = B;
   Filter3(ink).Fc = Fc(ink);  

%%	Filter resampled for fc = 125 Hz.
ink=3;
q=2; 	
FsNew=100*floor(Fs/q/100);
[B,A] = octdsgn(Fc(ink),FsNew,FiltOrd);
   Filter3(ink).FS = FsNew;
   Filter3(ink).A = A;
   Filter3(ink).B = B;
   Filter3(ink).Fc = Fc(ink);  

%%	Other filters should be fine
for ink=4:length(Fc)
   q=1;
   FsNew=Fs/q;
   [B,A] = octdsgn(Fc(ink),FsNew,FiltOrd);
   Filter3(ink).FS = FsNew;
   Filter3(ink).A = A;
   Filter3(ink).B = B;
   Filter3(ink).Fc = Fc(ink);   
end   
