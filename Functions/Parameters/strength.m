function G = strength(px, ps, fs, tDirect)
    % Calculates the "acoustic strength" from one signal px and its
    % anechoic version ps, measured at a distance of s = 10 m from the
    % sound source.
    % tDirect is the duration of direct sound in milliseconds.
    
    % TODO: implement fractional octave band filters
    
    if nargin < 4% if no tDirect is given...
        tDirect = 7;% ...assume a direct sound duration of 7 ms
    end
    
    % "Converts" tDirect to samples by multiplying the sampling frequency
    nDirect = round(tDirect / 1000 * fs);
    
    if length(ps) < nDirect
        error('Length of the anechoic reference signal ''ps'' is shorter than duration of the direct sound ''tDirect''.')
    else
        ps = ps(1:nDirect);
        
        % Assumes that an omnidirectional speaker was used so ? = 1
        G = 10 * log10( sum(px.^2) / sum(ps.^2) ) - 10 * log10( 4*pi*10.^2 );
    end
end
