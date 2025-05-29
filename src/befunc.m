function fx = befunc(x)
% Bi-elliptic radius objective function.
% Required by bielliptic.m

global mu ri rf

% Compute semi-major axes for transfer orbits.
sma1 = (ri + x) / 2;
sma2 = (x + rf) / 2;

% Compute velocities (km/s).
vi   = sqrt(mu / ri);
vt1a = sqrt((2*mu/ri) - (mu/sma1));
vt1b = sqrt((2*mu/x) - (mu/sma1));
vt2b = sqrt((2*mu/x) - (mu/sma2));
vt2c = sqrt((2*mu/rf) - (mu/sma2));
vf   = sqrt(mu / rf);

% Delta-v contributions.
dva = abs(vt1a - vi);
dvb = abs(vt2b - vt1b);
dvc = abs(vf - vt2c);

fx = dva + dvb + dvc;
end
