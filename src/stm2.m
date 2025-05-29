function [rf, vf, stm] = stm2(mu, tau, ri, vi)
% Compute the two-body state transition matrix using Shepperd's method.
% Inputs:
%   mu  - gravitational constant (km^3/s^2)
%   tau - propagation time interval (s)
%   ri, vi - initial ECI position and velocity vectors
% Outputs:
%   rf, vf - final ECI position and velocity vectors
%   stm - 6x6 state transition matrix

tol = 1.0e-8;
n0  = dot(ri, vi);
r0  = norm(ri);
beta = (2 * mu / r0) - dot(vi, vi);

u = 0;
if beta ~= 0
    umax = 1 / sqrt(abs(beta));
else
    umax = 1.0e24;
end
umin = -umax;

if beta > 0
    p = 2 * pi * mu * beta^(-1.5);
    n = fix((tau + 0.5 * p - 2*n0/beta) / p);
    delu = 2 * n * pi * beta^(-2.5);
else
    delu = 0;
end

tsav = 1.0e99;
niter = 0;

% Kepler iteration loop
while true
    niter = niter + 1;
    q = beta * u^2 / (1 + beta*u^2);
    u0 = 1 - 2*q;
    u1 = 2 * u * (1 - q);
    
    % Continued fraction iteration
    n_cf = 0; l = 3; d = 15; k = -9; a = 1; b = 1; g = 1;
    while true
        gsav = g;
        k = -k;
        l = l + 2;
        d = d + 4*l;
        n_cf = n_cf + (1 + k)*l;
        a = d / (d - n_cf * a * q);
        b = (a - 1) * b;
        g = g + b;
        if abs(g - gsav) < tol
            break;
        end
    end
    
    uu = (16/15) * u1^5 * g + delu;
    u2 = 2 * u1^2;
    u1 = 2 * u0 * u1;
    u0 = 2 * u0^2 - 1;
    u3 = beta * uu + (u1 * u2)/3;
    
    r1 = r0*u0 + n0*u1 + mu*u2;
    t  = r0*u1 + n0*u2 + mu*u3;
    dtdu = 4 * r1 * (1 - q);
    
    if abs(t - tsav) < tol, break; end
    usav = u; tsav = t;
    terr = tau - t;
    if abs(terr) < abs(tau)*tol, break; end
    du = terr / dtdu;
    
    if du < 0
        umax = u;
        u = u + du;
        if u < umin, u = 0.5*(umin + umax); end
    else
        umin = u;
        u = u + du;
        if u > umax, u = 0.5*(umin + umax); end
    end
    
    if abs(u - usav) < tol, break; end
    if niter > 20, break; end
end

% Compute Lagrange coefficients
fm  = -mu * u2 / r0;
ggm = -mu * u2 / r1;
f   = 1 + fm;
g   = r0*u1 + n0*u2;
ff  = -mu * u1 / (r0*r1);
gg  = 1 + ggm;

% Final state vectors
rf = f*ri + g*vi;
vf = ff*ri + gg*vi;

uu = g*u2 + 3*mu*uu;
a0 = mu / r0^3;
a1 = mu / r1^3;

% Preallocate stm and compute auxiliary matrix m.
stm = zeros(6,6);
m = zeros(3,3);
m(1,1) = ff * (u0/(r0*r1) + 1/r0^2 + 1/r1^2);
m(1,2) = (ff*u1 + (ggm/r1)) / r1;
m(1,3) = ggm*u1 / r1;
m(2,1) = -(ff*u1 + (fm/r0)) / r0;
m(2,2) = -ff*u2;
m(2,3) = -ggm*u2;
m(3,1) = fm*u1 / r0;
m(3,2) = fm*u2;
m(3,3) = g*u2;
m(1,1) = m(1,1) - a0*a1*uu;
m(1,3) = m(1,3) - a1*uu;
m(3,1) = m(3,1) - a0*uu;
m(3,3) = m(3,3) - uu;

% Fill state transition matrix
for i = 1:2
    xval = 2*i - 3;
    ib = 3*(2 - i);
    for j = 1:2
        jb = 3*(j - 1);
        for ii = 1:3
            t1 = rf(ii)*m(i,j) + vf(ii)*m(i+1,j);
            t2 = rf(ii)*m(i,j+1) + vf(ii)*m(i+1,j+1);
            for jj = 1:3
                stm(ii+ib, jj+jb) = xval * (t1*ri(jj) + t2*vi(jj));
            end
        end
    end
end

% Add Lagrange coefficients to stm submatrices.
for i = 1:3
    stm(i,i)     = stm(i,i)     + f;
    stm(i,i+3)   = stm(i,i+3)   + g;
    stm(i+3,i)   = stm(i+3,i)   + ff;
    stm(i+3,i+3) = stm(i+3,i+3) + gg;
end
end
