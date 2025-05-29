function [rf, vf] = twobody2(mu, tau, ri, vi)
% Solve the two-body initial value problem using Shepperd's method.
% Inputs:
%   mu  - gravitational constant (km^3/s^2)
%   tau - time interval (s)
%   ri, vi - initial position and velocity vectors
% Outputs:
%   rf, vf - final position and velocity vectors

tolerance = 1.0e-12;
u = 0;
imax = 20;
umax = realmax;
umin = -realmax;
orbits = 0;
tdesired = tau;
threshold = tolerance * abs(tdesired);
r0 = norm(ri);
n0 = dot(ri, vi);
beta = 2*(mu/r0) - dot(vi, vi);

if beta ~= 0
    umax = 1 / sqrt(abs(beta));
    umin = -1 / sqrt(abs(beta));
end

if beta > 0
    orbits = beta*tau - 2*n0;
    orbits = 1 + (orbits*sqrt(beta))/(pi*mu);
    orbits = floor(orbits/2);
end

for i = 1:imax
    q = beta * u^2 / (1 + beta*u^2);
    r_cf = 1; n_cf = 0; l = 1; s = 1; d = 3; k = -5;
    gcf = 1; gold = 0;
    while gcf ~= gold
        k = -k;
        l = l + 2;
        d = d + 4*l;
        n_cf = n_cf + (1+k)*l;
        r_cf = d / (d - n_cf * r_cf * q);
        s = (r_cf - 1) * s;
        gold = gcf;
        gcf = gold + s;
    end
    
    h0 = 1 - 2*q;
    h1 = 2*u*(1 - q);
    u0 = 2*h0^2 - 1;
    u1 = 2*h0*h1;
    u2 = 2*h1^2;
    u3 = 2*h1*u2*gcf/3;
    
    if orbits ~= 0
        u3 = u3 + 2*pi*orbits/(beta*sqrt(beta));
    end
    
    r1 = r0*u0 + n0*u1 + mu*u2;
    dt = r0*u1 + n0*u2 + mu*u3;
    slope = 4*r1 / (1 + beta*u^2);
    terror = tdesired - dt;
    
    if abs(terror) < threshold
        break;
    end
    
    if (i > 1) && (u == uold || dt == dtold)
        break;
    end
    uold = u; dtold = dt;
    ustep = terror / slope;
    if ustep > 0
        umin = u;
        u = u + ustep;
        if u > umax, u = (umin+umax)/2; end
    else
        umax = u;
        u = u + ustep;
        if u < umin, u = (umin+umax)/2; end
    end
    if i == imax
        fprintf('\nMaximum iterations reached in twobody2.');
    end
end

f = 1 - (mu/r0)*u2;
gg = 1 - (mu/r1)*u2;
g = r0*u1 + n0*u2;
ff = -mu*u1/(r0*r1);

% Vectorized final state computation.
rf = f * ri + g * vi;
vf = ff * ri + gg * vi;
end
