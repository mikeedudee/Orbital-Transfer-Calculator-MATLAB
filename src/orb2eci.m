function [r, v] = orb2eci(mu, oev)
% Convert classical orbital elements to ECI state vectors.
% Inputs:
%   mu  - gravitational constant (km^3/s^2)
%   oev - [sma, eccentricity, inclination, arg perigee, RAAN, true anomaly]
% Outputs:
%   r - position vector (km)
%   v - velocity vector (km/s)

sma    = oev(1);
ecc    = oev(2);
inc    = oev(3);
argper = oev(4);
raan   = oev(5);
tanom  = oev(6);

slr = sma * (1 - ecc^2);
rm  = slr / (1 + ecc*cos(tanom));
arglat = argper + tanom;

% Precompute trigonometric terms.
sarglat = sin(arglat);
carglat = cos(arglat);
c4 = sqrt(mu / slr);
c5 = ecc*cos(argper) + carglat;
c6 = ecc*sin(argper) + sarglat;
sinc = sin(inc);
cinc = cos(inc);
sraan = sin(raan);
craan = cos(raan);

% Position vector in ECI.
r = [rm * (craan * carglat - sraan * cinc * sarglat);
     rm * (sraan * carglat + cinc * sarglat * craan);
     rm * sinc * sarglat];

% Velocity vector in ECI.
v = [-c4 * (craan * c6 + sraan * cinc * c5);
     -c4 * (sraan * c6 - craan * cinc * c5);
      c4 * c5 * sinc];
end
