function pviniz(tof, r1, v1, dv1, dv2)
% Primer vector initialization.
% Required by primer.m

global mu pvi pvdi

% Compute normalized dv vectors.
pvi = dv1 / norm(dv1);
pvf = dv2 / norm(dv2);

% Compute state transition matrix at time tof.
[~, ~, stm] = stm2(mu, tof, r1, v1);

% Extract submatrices.
stm11 = stm(1:3,1:3);
stm12 = stm(1:3,4:6);

% Compute initial primer derivative vector using matrix division.
pvdi = stm12 \ (pvf' - stm11 * pvi');
end
