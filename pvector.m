function [pvm, pvdm] = pvector(ri, vi, x)
% Compute the primer vector and its derivative magnitude.
% Required by primer.m

global mu pvi pvdi

% Compute state transition matrix at time x.
[~, ~, stm] = stm2(mu, x, ri, vi);

% Evaluate the primer vector fundamental equation.
ppdot = stm * [pvi'; pvdi];

% Extract primer vector and its derivative.
pv  = ppdot(1:3);
pvd = ppdot(4:6);

% Compute magnitudes.
pvm  = norm(pv);
pvdm = dot(pv, pvd) / pvm;
end
