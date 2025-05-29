%% Bi-Elliptic Orbit Transfer Analysis and 3D Visualization
clear; clc;
global mu req ri rf pvi pvdi dtr

% Set conversion factors and constants
om_constants;  % loads mu, req, etc.

    global alti altf   % Declare globals

    global dtr rtd mu mmu smu omega req flat j2 aunit
    om_constants;
    format long;
    EarthRadius           = ((2*6378.1370)+6356.7523)/3;
    GravitationalConstant = 6.67430e-11;
    EarthMass             = (5.9722+0.0006)*1e24;
    u                     = GravitationalConstant * EarthMass;
    dtr   = pi / 180.0;

% Request inputs
fprintf('<strong>One Tangent Transfer\n</strong>');

disp(['Enter the initial altitude [' 8 '(km)]' 8 ':'])
alti = get_positive_input('> ');

disp(['Enter the final altitude [' 8 '(km)]' 8 ':'])
altf = get_positive_input('> ');

% Use consistent naming if you want to get the transfer ellipse altitude from the user:
fprintf('Enter transfer ellipse semi-major axis altitude (km):\n');
transferAlt = get_between_altitudes('> ');   % new name for clarity


fprintf('Set the orbits incliniation:\n');
incl = get_nomorethan('> ');

fprintf('Set the orbits right ascension of the ascending node:\n');
raan = input('> ');

    % Use user inputs in calculations:
    orbitI = (EarthRadius + alti) * 1000;
    orbitF = (EarthRadius + altf) * 1000;
    Atx    = (transferAlt + EarthRadius) * 1000;

    e = 1 - (orbitI/Atx);
    v = acosd((((Atx * (1 - e^2))/orbitF)-1)/e);
    a = atand((e * sind(v)) / (1 + e * cosd(v)));
    Vi1 = sqrt(u/orbitI);
    Vf2 = sqrt(u/orbitF);
    Vtx1 = sqrt(u * ((2/orbitI)-(1/Atx)));
    Vtx2 = sqrt(u * ((2/orbitF)-(1/Atx)));
    DV1  = abs(Vtx1 - Vi1);
    DV2  = abs(sqrt(Vtx2^2 + Vf2^2 - 2 * Vtx2 * Vf2 * cosd(a)));
    DVT = DV1 + DV2;
    E = acos(((orbitF/Atx)-1)/-e);  % E in radians
    TOF = (E - (e * sin(E))) * sqrt(Atx^3/u);
    TOFs = TOF/3600;
    period1 = sqrt((4*pi^2*orbitI^3)/u)/3600;
    period2 = sqrt((4*pi^2*orbitF^3)/u)/3600;

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
print_results(period1, period2, orbitI, orbitF, e, v, a, Vi1, Vf2, Vtx1, Vtx2, DV1, DV2, DVT, E, TOF, Atx);

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
deltat3  = period3 / 300;
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
rx = [rp1_x, rp2_x, rp3_x];
ry = [rp1_y, rp2_y, rp3_y];
rz = [rp1_z, rp2_z, rp3_z];

% Optionally, plot the full path for reference
%plot3(rx, ry, rz, 'k--', 'LineWidth', 1.5);

% Initialize a marker for the rocket (red circle)
rocketHandle = plot3(rx(1), ry(1), rz(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
secondObjectHandle = plot3(rx(1), ry(1), rz(1), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'c');
meeting_index = 0;
offset = meeting_index - 1;

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
legend('Rocket','Target Object','Earth','Final Orbit','Initial Orbit','Transfer Ellipse', 'Location','best');
title('One Tangent Transfer', 'FontSize', 16);
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

nRocket = length(rx);    % Total points in rocket's trajectory
nOrbit  = length(rp3_x);   % Total points for orbit data
i = 1;
while true
    % For the rocket:
    if i <= nRocket
        % Use the rocket's own trajectory until its data is exhausted.
        set(rocketHandle, 'XData', rx(i), 'YData', ry(i), 'ZData', rz(i));
    else
        % Once the rocket reaches its max index, start using the orbit data.
        % Wrap the index to cycle continuously through rp3 data.
        idx_r = mod(i - nRocket - 1, nOrbit) + 1;
        set(rocketHandle, 'XData', rp3_x(idx2), 'YData', rp3_y(idx2), 'ZData', rp3_z(idx2));
    end
    
    % For the second object, always use the orbit (rp3) data with wrapping.
    idx2 = mod(i - 1 + offset, nOrbit) + 1;
    set(secondObjectHandle, 'XData', rp3_x(idx2), 'YData', rp3_y(idx2), 'ZData', rp3_z(idx2));
    
    drawnow;
    %pause(0.05); % Adjust pause for desired animation speed
    
    i = i + 1;
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

function val = get_between_altitudes(prompt)
    global alti altf
    while true
        val = input(prompt);
        if isnumeric(val) && (val < altf) && (val > alti)
            break;
        else
            fprintf('Value must be between the initial and final altitudes.\n');
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

function print_results(period1, period2, orbitI, orbitF, e, v, a, Vi1, Vf2, Vtx1, Vtx2, DV1, DV2, DVT, E, TOF, Atx);
    fprintf('<strong>One Tangent Analysis Results\n</strong>');
    fprintf('--------------------------------------------\n');
    fprintf('Initial orbit altitude:      %14.8f m\n', orbitI);
    fprintf('Final orbit altitude:        %14.8f m\n', orbitF);
    fprintf('Initial velocity:            %14.8f m/s\n', Vi1);
    fprintf('Final velocity:              %14.8f m/s\n', Vf2);
    fprintf('Orbital Period 1st Orbit:    %14.8f Hr\n', period1);
    fprintf('Orbital Period 2nd Orbit:    %14.8f Hr\n', period2);
    fprintf('Transfer ellipse (alt SMA):  %14.8f km\n', Atx);
    fprintf('Eccentricity:                %14.8f \n', e);
    fprintf('True Anomaly @ 2nd burn:     %14.8f degrees\n', v);
    fprintf('Path angle @ 2nd burn:       %14.8f degrees\n', a);
    fprintf('Eccentricity Anomaly:        %14.8f \n', E);
    fprintf('Time-of-Flight E:            %14.8f Hr\n', TOF);
    fprintf('Velocity transfer Initial:   %14.8f m/s\n', Vtx1);
    fprintf('Velocity transfer Final:     %14.8f m/s\n', Vtx2);
    fprintf('Delta V Initial:             %14.8f m/s\n', DV1);
    fprintf('Delta V Final:               %14.8f m/s\n', DV2);
    fprintf('Delta V Total:               %14.8f m/s\n', DVT);
end