function om_constants
% Define astrodynamic and utility constants.

global dtr rtd mu mmu smu omega req flat j2 aunit

dtr   = pi / 180.0;
rtd   = 180.0 / pi;
mu    = 398600.436233;               % Earth gravitational constant (km^3/s^2)
mmu   = 4902.800076;                 % Moon gravitational constant (km^3/s^2)
smu   = 132712440040.944;            % Sun gravitational constant (km^3/s^2)
omega = 7.292115486e-5;              % Earth inertial rotation rate (rad/s)
req   = ((2*6378.1370)+6356.7523)/3; % Earth equatorial radius (km)
flat  = 1.0 / 298.257;               % Earth flattening factor
j2    = 0.00108263;                  % Earth oblateness coefficient
aunit = 149597870.691;               % Astronomical unit (km)
end
