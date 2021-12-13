function y = magfit_setspacing(x,o)
    % Constants for GRIN lens
    % x = optical axis distance
    % a = grin to relay spacing
    % o = offset
    
    z = x+o
    
    l = 2.46;
    n0 = 1.635;
    P = 0.5;
    
    a = 19.0516;
    
    
    di0 = ImageDistanceGRIN(n0, l, P, z);
    dGRIN = 13.4 + a; % convert to distance to relay lens plane
                        %13.4 mm from the doublet tube front to the lens plane
    grinMag = -di0 ./ z;
    do1 = dGRIN - di0;
    di1 = ((1/30) - (1 ./ do1)) .^ -1;
    relay1Mag = -di1 ./ do1;
    do2 = 16.6 - di1; % 16.6 mm is distance between principle planes of relay lenses
    di2 = ((1/40) - (1 ./ do2)) .^ -1;
    relay2Mag = -di2 ./ do2;
    y = grinMag .* relay1Mag .* relay2Mag;

end