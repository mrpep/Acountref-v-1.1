function G = strength(px, ps, fs, tDirect, gamma)
    % Calculates the "acoustic strength" from one signal px and its
    % anechoic version ps, measured at a distance of s = 10 m from the
    % sound source.
    % tDirect is the duration of direct sound in milliseconds.
    
    if nargin < 5% if no gamma is given...
        gamma = 1;% ... assume that an omnidirectional loudspeaker was used for the measurement.
        if nargin < 4% if no tDirect is given either...
            tDirect = 7;% ...assume a direct sound duration of 7 ms
        end
    end
    
    % "Converts" tDirect to samples by multiplying the sampling frequency
    nDirect = round(tDirect / 1000 * fs);
    
    if length(ps) < nDirect
        error('Length of the anechoic reference signal ''ps'' is shorter than duration of the direct sound ''tDirect''.')
    else
        ps = ps(1:nDirect);
        
        % Assumes that an omnidirectional speaker was used so gamma = 1
        G = 10 * log10( trapz(px.^2) / ( gamma * trapz(ps.^2) ) ) - 10 * log10( 4*pi*10.^2 );
    end
end
