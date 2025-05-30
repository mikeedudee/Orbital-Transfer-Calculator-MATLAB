clear; clc;
efficiency = 1;
global Earth_gravitational_constant Earth_radius degrees_to_radians radians_to_degrees altf;

cvc = contantsvalues_convertions(); % load the constants of the function

% Request input from the user
fprintf('<strong>Bi-Elliptic Orbit Transfer\n</strong>');

disp(['Enter the initial altitude [' 8 '(km)]' 8 ':'])
alti = get_positive_input('> ');

disp(['Enter the final altitude [' 8 '(km)]' 8 ':'])
altf = get_positive_input('> ');

disp(['Enter the intermediate apoapsis altitude [' 8 '(km)]' 8 ':'])
alt_intermediate = get_nolessthan('> ');

fprintf('Set the orbits inclination:\n');
incl = get_nomorethan('> ');

disp(['Right Ascension of the Ascending Node (RAAN) in degrees [' 8 ']' 8 ':'])
raan = input('> '); %Right Ascension of the Ascending Node (RAAN) in degrees

%% Getting the computations
[radius_initial, radius_final, velocity_initial_orbit, velocity_final_orbit, ...
    semi_major_axis_first, eccentricity_first, periapsis_velocity_i, apoapsis_velocity_i, ...
    semi_major_axis_second, eccentricity_second, periapsis_velocity_ii, apoapsis_velocity_ii, ...
    time_initial_orbit, time_intermediate_orbit, time_final_orbit, ...
    total_orbit_period, transfer_time_first, transfer_time_second, total_transfer_time, ...
    deltaV_1, deltaV_2,  deltaV_3, deltaV_total, efficiency, radius_intermediate] = bielliptic_orbital_parameters(alti, altf, alt_intermediate, incl, efficiency);

%% Display the results
% Precompute converted units
radius_initial_km       = radius_initial / 1e3;
radius_intermediate_km  = radius_intermediate / 1e3; % Convert intermediate altitude to km
radius_final_km         = radius_final / 1e3;
semimajor_axis_first_km = semi_major_axis_first / 1e3;
semimajor_axis_secon_km = semi_major_axis_second / 1e3;
time_initial_hr         = time_initial_orbit / 3600;
time_intermediate_hr    = time_intermediate_orbit / 3600;
time_final_hr           = time_final_orbit / 3600;
total_period_hr         = total_orbit_period / 3600;
transfer_time_first_hr  = transfer_time_first / 3600;
transfer_time_second_hr = transfer_time_second / 3600;
transfer_time           = total_transfer_time / 3600;
raan_deg                = raan * cvc.radians_to_degrees;

% Labels and values
labels = {
    'Initial Orbital Radius',       radius_initial_km,       'km';
    'Intermediate Orbital Radius',  radius_intermediate_km,  'km';
    'Final Orbital Radius',         radius_final_km,         'km';
    'Initial Orbital Velocity',     velocity_initial_orbit,  'm/s';
    '1st Transfer Semi-major Axis', semimajor_axis_first_km, 'km';
    '1st Transfer Eccentricity',    eccentricity_first,      '';
    '1st Periapsis Velocity',       periapsis_velocity_i,    'm/s';
    '1st Apoapsis Velocity',        apoapsis_velocity_i,     'm/s';
    '2nd Transfer Semi-major Axis', semimajor_axis_secon_km, 'km';
    '2nd Transfer Eccentricity',    eccentricity_second,     '';
    '2nd Periapsis Velocity',       periapsis_velocity_ii,   'm/s';
    '2nd Apoapsis Velocity',        apoapsis_velocity_ii,    'm/s';
    'Final Orbit Velocity',         velocity_final_orbit,    'm/s';
    'Initial Orbit Time',           time_initial_hr,         'hour(s)';
    'Intermediate Orbit Time',      time_intermediate_hr,    'hour(s)';
    'Final Orbit Time',             time_final_hr,           'hour(s)';
    'Total Orbit Period',           total_period_hr,         'hour(s)';
    '1st Transfer Time',            transfer_time_first,     'hour(s)';
    '2nd Transfer Time',            transfer_time_second,    'hour(s)';
    'Total Transfer Time',          total_transfer_time,     'hour(s)';
    'Delta V for First Transfer',   deltaV_1,                'm/s';
    'Delta V for Second Transfer',  deltaV_2,                'm/s';
    'Delta V for Third Transfer',   deltaV_3,                'm/s';
    'Total Delta V',                deltaV_total,            'm/s';
};


% Display using a loop
fprintf('\n--------------------------------------------\n');
fprintf('<strong>Hohmann Transfer Analysis Results\n</strong>');
fprintf('--------------------------------------------\n');
    for k = 1:size(labels,1)
        fprintf('%s: %.4f %s\n', labels{k,1}, labels{k,2}, labels{k,3});
    end

    if efficiency == 0
        fprintf(2,'\nThe Hohmann method is efficient for this kind of transfer,\n');
        fprintf(2,'use it instead.\n');
    end

% Generate the figure
figure_main()
title('BI-ELLIPTIC TRANSFER', 'FontSize', 15);

% Draw the Earth sphere
earth_sphere(); 

% Conversion of inclination and RAAN from degrees to radians
i_rad = incl * cvc.degrees_to_radians;

% Set the argument of periapsis to zero for Hohmann transfer
omega = 0; % [rad]

% True anomaly range for full orbit
theta = linspace(0, 2*pi, 1000); 

% Preallocate position array
r_eci = zeros(3, length(theta));

% Compute the orbits
[initial_orbit, first_half, second_half, final_orbit] = orbits(raan, i_rad, omega, radius_initial, ...
radius_final, theta, incl, semi_major_axis_first, eccentricity_first, semi_major_axis_second, eccentricity_second, true);

%compute_to_eci();

Earth_radius_m = cvc.Earth_radius * 1e3;

% Call the plotting function to draw the orbits and the rocket marker
% Initialize the rocket marker at the first point of the initial orbit
[rocket_marker] = plotting(Earth_radius_m, initial_orbit, first_half, second_half, final_orbit, true);

% Call the animation function
animation(first_half, second_half, final_orbit, initial_orbit, Earth_radius_m, rocket_marker, true);

%% Local Function for validating the user input
function val = get_positive_input(prompt)
    while true
        val = input(prompt);
        if isnumeric(val) && (val > 0) && isscalar(val)
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

function val = get_nolessthan(prompt)
    global altf;
    while true
        val = input(prompt);
        if isnumeric(val) && (val > altf)
            break;
        else
            fprintf('Value must not be less than the target orbit.\n');
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