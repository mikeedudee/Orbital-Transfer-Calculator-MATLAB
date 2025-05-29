%% Bi-Elliptic Orbit Transfer Analysis and 3D Visualization
clear; clc;
global mu req ri rf pvi pvdi dtr

% Set conversion factors and constants
om_constants;  % loads mu, req, etc.

% Request inputs
fprintf('<strong>Bi-Elliptic Orbit Transfer\n</strong>');

disp(['Enter the initial altitude [' 8 '(km)]' 8 ':'])
alti = get_positive_input('> ');

disp(['Enter the final altitude [' 8 '(km)]' 8 ':'])
altf = get_positive_input('> ');

fprintf('Set the orbits incliniation:\n');
incl = get_nomorethan('> ');

fprintf('Set the orbits right ascension of the ascending node:\n');
raan = input('> ');
        % initial and final orbit radii (kilometers)
        ri = alti + req;
        rf = altf + req;
        % compute optimum bi-elliptic radius
        xmin              = rf;
        xmax              = 100.0 * rf;
        options           = optimset('TolFun', 1.0e-6, 'TolX', 1.0e-6);
        [x, fx, exitflag] = fminbnd('befunc', xmin, xmax, options);

        % intermediate orbit radius and altitude (kilometers)
        rb   = x;
        altb = rb - req;

% semimajor axes of the elliptical transfer orbits (kilometers)
smat1 = (ri + rb) / 2.0;
smat2 = (rb + rf) / 2.0;

% compute first ellipse radii (kilometers) and eccentricity (non-dimensional)
rpt1  = req + alti;
rat1  = req + altb;
ecct1 = (rat1 - rpt1) / (rat1 + rpt1);

% compute second ellipse radii (kilometers) and eccentricity (non-dimensional)
rpt2  = req + altf;
rat2  = req + altb;
ecct2 = (rat2 - rpt2) / (rat2 + rpt2);
vi    = sqrt(mu / ri);                        % initial circular orbit velocity (kilometers/second)
vt1a  = sqrt((2.0 * mu / ri) - (mu / smat1)); % first ellipse perigee velocity (kilometers/second)
vt1b  = sqrt((2.0 * mu / rb) - (mu / smat1)); % first ellipse apogee velocity (kilometers/second)
vt2b  = sqrt((2.0 * mu / rb) - (mu / smat2)); % second ellipse apogee velocity (kilometers/second)
vt2c  = sqrt((2.0 * mu / rf) - (mu / smat2)); % second ellipse perigee velocity (kilometers/second)
vf    = sqrt(mu / rf);                        % final circular orbit velocity (kilometers/second)

% compute delta-v contributions (kilometers/second)
dva = abs(vt1a - vi);
dvb = abs(vt2b - vt1b);
dvc = abs(vf - vt2c);

% compute individual and total transfer time (seconds)
tof1 = pi * sqrt(smat1^3 / mu);
tof2 = pi * sqrt(smat2^3 / mu);
tof  = tof1 + tof2;

% print results
clc; home;

% Display results
print_results(alti, altf, altb, vi, vt1a, vt1b, vt2b, vt2c, vf, ecct1, ecct2, ...
    dva, dvb, dvc, tof1, tof2, tof);

% Set desired inclination and RAAN in radians
inclination = incl * dtr;   % 30 degrees converted to radians
raan        = raan * dtr;   % 40 degrees converted to radians

% Update the initial orbit elements: [a, e, i, ω, Ω, θ]
oevi(1)  = ri;           % semimajor axis (or orbit radius for circular orbit)
oevi(2)  = 0.0;           % eccentricity (circular orbit)
oevi(3)  = inclination;   % inclination
oevi(4)  = 0.0;           % argument of perigee (can be set as needed)
oevi(5)  = raan;          % RAAN
oevi(6)  = 0.0;           % true anomaly
[ri, vi] = orb2eci(mu, oevi);

% orbital elements and state vector at perigee of the first ellipse
oevti1(1) = smat1;
oevti1(2) = ecct1;
oevti1(3) = inclination;
oevti1(4) = 0.0;
oevti1(5) = raan;
oevti1(6) = 0.0;
[rtp1, vtp1] = orb2eci(mu, oevti1);

% orbital elements and state vector at apogee of the first ellipse
oevtf2(1) = smat1;
oevtf2(2) = ecct1;
oevtf2(3) = inclination;
oevtf2(4) = 0.0;
oevtf2(5) = raan;
oevtf2(6) = 180.0 * dtr;
[rta1, vta1] = orb2eci(mu, oevtf2);

% orbital elements and state vector at apogee of the second ellipse
oevti2(1) = smat2;
oevti2(2) = ecct2;
oevti2(3) = inclination;
oevti2(4) = 0.0;
oevti2(5) = raan;
oevti2(6) = 180.0 * dtr;
[rta2, vta2] = orb2eci(mu, oevti2);

% orbital elements and state vector at perigee of the second ellipse
oevtf2(1) = smat2;
oevtf2(2) = ecct2;
oevtf2(3) = inclination;
oevtf2(4) = 0.0;
oevtf2(5) = 0.0;
oevtf2(6) = 0.0;
[rtp2, vtp2] = orb2eci(mu, oevtf2);

% orbital elements and state vector of final circular orbit
oevf(1) = rf;
oevf(2) = 0.0;
oevf(3) = inclination;
oevf(4) = 0.0;
oevf(5) = raan;
oevf(6) = 0.0;
[rf, vf] = orb2eci(mu, oevf);

% compute orbital periods (seconds)
period1  = 2.0 * pi * oevi(1) * sqrt(oevi(1) / mu);
period2  = 2.0 * pi * oevti1(1) * sqrt(oevti1(1) / mu);
period3  = 2.0 * pi * oevti2(1) * sqrt(oevti2(1) / mu);
period4  = 2.0 * pi * oevf(1) * sqrt(oevf(1) / mu);
deltat1  = period1 / 300;
simtime1 = -deltat1;
deltat2  = 0.5 * period2 / 300;
simtime2 = -deltat2;
deltat3  = 0.5 * period3 / 300;
simtime3 = -deltat3;
deltat4  = period4 / 300;
simtime4 = -deltat4;

for i = 1:1:301
    
    simtime1 = simtime1 + deltat1;
    simtime2 = simtime2 + deltat2;
    simtime3 = simtime3 + deltat3;
    simtime4 = simtime4 + deltat4;
    
    % initial orbit position vector (ER)
    [rwrk, vwrk] = twobody2 (mu, simtime1, ri, vi);
    rp1_x(i) = rwrk(1) / req;
    rp1_y(i) = rwrk(2) / req;
    rp1_z(i) = rwrk(3) / req;
    
    % first transfer orbit position vector (ER)
    [rwrk, vwrk] = twobody2 (mu, simtime2, rtp1, vtp1);
    rp2_x(i) = rwrk(1) / req;
    rp2_y(i) = rwrk(2) / req;
    rp2_z(i) = rwrk(3) / req;
    
    % second transfer orbit position vector (ER)
    [rwrk, vwrk] = twobody2 (mu, simtime3, rta2, vta2);
    
    rp3_x(i) = rwrk(1) / req;
    rp3_y(i) = rwrk(2) / req;
    rp3_z(i) = rwrk(3) / req;    
    
    % final orbit position vector (ER)
    [rwrk, vwrk] = twobody2 (mu, simtime4, rf, vf);
    rp4_x(i) = rwrk(1) / req;
    rp4_y(i) = rwrk(2) / req;
    rp4_z(i) = rwrk(3) / req;    
    
end

% create axes vectors
xaxisx = [1 1.5];
xaxisy = [0 0];
xaxisz = [0 0];

yaxisx = [0 0];
yaxisy = [1 1.5];
yaxisz = [0 0];

zaxisx = [0 0];
zaxisy = [0 0];
zaxisz = [1 1.5];

figure (1);
hold on;
grid on;

% Combine the trajectory points from all phases into one continuous set
rx = [rp1_x, rp2_x, rp3_x, rp4_x];
ry = [rp1_y, rp2_y, rp3_y, rp4_y];
rz = [rp1_z, rp2_z, rp3_z, rp4_z];

% Optionally, plot the full path for reference
%plot3(rx, ry, rz, 'k--', 'LineWidth', 1.5);

% Initialize a marker for the rocket (red circle)
rocketHandle = plot3(rx(1), ry(1), rz(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'y');

% plot earth
[x y z] = sphere(24);
h       = surf(x, y, z);
map     = [0 0 0.3
           0 0 0.4
           0 0 0.5
           0 0 0.6
           0 0 0.8
           0 0 1.0];
colormap(map);
set (h, 'edgecolor', [1 1 1]);
shading interp;  % Smooth the surface appearance

% plot coordinate system axes
plot3(xaxisx, xaxisy, xaxisz, '-g', 'LineWidth', 1);
plot3(yaxisx, yaxisy, yaxisz, '-r', 'LineWidth', 1);
plot3(zaxisx, zaxisy, zaxisz, '-b', 'LineWidth', 1);

% plot initial orbit
plot3(rp1_x, rp1_y, rp1_z, '-r', 'LineWidth', 1.5);
plot3(rp1_x(1), rp1_y(1), rp1_z(1), 'ob');

% plot first transfer orbit
plot3(rp2_x, rp2_y, rp2_z, '-b', 'LineWidth', 1.5);
plot3(rp2_x(end), rp2_y(end), rp2_z(end), 'ob');

% plot second transfer orbit
%plot3(rp3_x, rp3_y, rp3_z, '-b', 'LineWidth', 1.5);
plot3(rp3_x(end), rp3_y(end), rp3_z(end), 'ob');

% plot final orbit
plot3(rp4_x, rp4_y, rp4_z, '-g', 'LineWidth', 1.5);
xlabel('X coordinate (ER)', 'FontSize', 12);
ylabel('Y coordinate (ER)', 'FontSize', 12);
zlabel('Z coordinate (ER)', 'FontSize', 12);
legend('Rocket','Earth','Final Orbit','Initial Orbit','Transfer Ellipse','Location','best');
title('Bi-elliptic Transfer: Initial, Transfer and Final Orbits', 'FontSize', 16);
axis equal;
view(50, 40);
rotate3d on;
print -depsc -tiff -r300 bielliptic1.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process first transfer ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dvi = (vtp1 - vi)';                 % initial deltav
dvf = (vta2 - vta1)';               % final deltav
npts = 300;                         % number of graphic data points
dt = tof1 / npts;                   % create primer vector and derivative data
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process second transfer ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dvi = (vta2 - vta1)';               % initial deltav
dvf = (vf - vtp2)';                 % final deltav
pviniz(tof2, rta2, vta2, dvi, dvf); % perform primer vector initialization
npts = 300;                         % number of graphic data points
dt = tof2 / npts;                   % create primer vector and derivative data

% Animate the rocket along the path
for i = 1:length(rx)
    set(rocketHandle, 'XData', rx(i), 'YData', ry(i), 'ZData', rz(i));
    drawnow;             % Update the figure
    %pause(0.05);        % Adjust pause time for desired animation speed
end

%% Local Function Definitions
function val = get_positive_input(prompt)
    while true
        val = input(prompt);
        if isnumeric(val) && (val > 0)
            break;
        else
            fprintf('Value must be positive.\n');
        end
    end
end

function val = get_nomorethan(prompt)
    while true
        val = input(prompt);
        if isnumeric(val) && (val < 360)
            break;
        else
            fprintf('Value must not be more than 360.\n');
        end
    end
end

function val = get_selection_input(prompt, valid_options)
    while true
        val = input(prompt);
        if ismember(val, valid_options)
            break;
        else
            fprintf('Invalid selection. Try again.\n');
        end
    end
end

function [a, ecc] = compute_ellipse(rp, ra)
    a = (rp + ra) / 2;
    ecc = (ra - rp) / (ra + rp);
end

function [dva, dvb, dvc] = compute_deltav(vi, vt1a, vt1b, vt2b, vt2c, vf)
    dva = abs(vt1a - vi);
    dvb = abs(vt2b - vt1b);
    dvc = abs(vf - vt2c);
end

function print_results(alti, altf, altb, vi, vt1a, vt1b, vt2b, vt2c, vf, ecct1, ecct2, dva, dvb, dvc, tof1, tof2, tof)
    fprintf('<strong>Hooman Transfer Analysis Results\n</strong>');
    fprintf('--------------------------------------------\n');
    fprintf('Initial orbit altitude:      %14.4f km\n', alti);
    fprintf('Final orbit altitude:        %14.4f km\n', altf);
    fprintf('Initial velocity:            %14.4f m/s\n', vi * 1000);
    fprintf('Transfer orbit velocity:     %14.4f m/s\n', vt1a * 1000);
    fprintf('Final velocity:              %14.4f m/s\n', vf * 1000);
    fprintf('Eccentricity (transfer):     %14.8f\n', ecct1);
    fprintf('Delta-v:                     %14.4f, %14.4f m/s\n', dva * 1000, dvb * 1000);
    fprintf('Total delta-v:               %14.4f m/s\n', (dva + dvb) * 1000);
    fprintf('Transfer time:               %14.4f hours\n', tof / 3600);
end
