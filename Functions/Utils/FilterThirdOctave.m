function [y,Fs,Fc] = FilterThirdOctave(sig, fe, B)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [y,Fc] = FilterThirdOctave(sig,fe)
%
% FUNCTION
%       calculation of 1/3 octave levels in the 29 bands
%       Band 1: fc = 25Hz, Band 29: fc = 16000Hz
%       or octave level in the 10 bands
%       Band 1: fc = 31.5Hz, Band 10: fc = 16000Hz
%
% INPUT
%       sig = signal (Pa)
%       fe = sampling frequency (Hz)
%       B = octave band (1) or 1/3 octave band (3)
% 
% OUTPUT
%       y =  Filter signal
%       Fc = central frequencies of third octave bands (Hz)
%       Fs = New Sampling frequency of each filter (Hz)       
%
% CREDITS
%   Third octave filters used in this function are coming from a program
%   from Christophe Couvreur which can be found in MatlabCentral (http://www.mathworks.com/matlabcentral/)
%   and modified by Aaron Hastings, Herrick Labs, from Purdue University 
%   in order to have better precision at low frequencies (subsampling)
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - 2009 %
%   Modified by Florent Masson - 2015  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if B==3 

% get IIR filters
Nbandes = 29;
StructFilt = ThirdOctave_Filters(fe); % structure of IIR for each band

Fc = zeros(1, Nbandes);
for k=1:Nbandes,
    Fc(k) = StructFilt(k).Fc;
    Fs(k) = StructFilt(k).FS;
end;
     
      
%% calculation per band
for i = 1:Nbandes,
    % resampling if necessary
    if fe ~= StructFilt(i).FS,
        sig_tmp = resample(sig, StructFilt(i).FS, fe);
    else
        sig_tmp = sig;
    end;
    y{i} = zeros(length(sig_tmp),1);
    % filtering
    y{i} = filter(StructFilt(i).B, StructFilt(i).A, sig_tmp);
end;
 
else
    
Nbandes = 9;
StructFilt = Octave_Filters(fe); % structure of IIR for each band

Fc = zeros(1, Nbandes);
for k=1:Nbandes,
    Fc(k) = StructFilt(k).Fc;
    Fs(k) = StructFilt(k).FS;
end;
     
      
%% calculation per band
for i = 1:Nbandes,
    % resampling if necessary
    if fe ~= StructFilt(i).FS,
        sig_tmp = resample(sig, StructFilt(i).FS, fe);
    else
        sig_tmp = sig;
    end;
    y{i} = zeros(length(sig_tmp),1);
    % filtering
    y{i} = filter(StructFilt(i).B, StructFilt(i).A, sig_tmp);
end;    
    
end


