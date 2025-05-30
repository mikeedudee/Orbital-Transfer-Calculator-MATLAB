clear; clc;
efficiency = 0;
global Earth_gravitational_constant Earth_radius degrees_to_radians radians_to_degrees;

cvc = contantsvalues_convertions(); % load the constants of the function

% Request input from the user
fprintf('<strong>One-Tangent Orbit Transfer\n</strong>');

disp(['Enter the initial altitude [' 8 '(km)]' 8 ':'])
alti = get_positive_input('> ');

disp(['Enter the target altitude [' 8 '(km)]' 8 ':'])
altf = get_positive_input('> ');

disp(['Right Ascension of the Ascending Node (RAAN) in degrees [' 8 ']' 8 ':'])
raan = input('> '); %Right Ascension of the Ascending Node (RAAN) in degrees

%% Getting the computations
[radius_initial, radius_final, velocity_initial_orbit, velocity_final_orbit, ...
    semimajor_axis, periapsis_velocity, apoapsis_velocity, ...
    time_final_orbit, time_initial_orbit, total_orbit_period, ...
    deltaV_1, deltaV_2, deltaV_total, transfer_time, efficiency, Phase_AngleD] = onetangent_orbital_parameters(alti, altf, efficiency);

%% Display the results
% Precompute converted units
radius_initial_km     = radius_initial / 1e3;
radius_final_km       = radius_final / 1e3;
semimajor_axis_km     = semimajor_axis / 1e3;
time_initial_hr       = time_initial_orbit / 3600;
time_final_hr         = time_final_orbit / 3600;
total_period_hr       = total_orbit_period / 3600;
transfer_time_hr      = transfer_time;
raan_deg              = raan * cvc.radians_to_degrees;

% Labels and values (dictionary-like structure)
labels = {
    'Earth Radius',               cvc.Earth_radius,          'km';
    'Initial Orbital Radius',     radius_initial_km,         'km';
    'Final Orbital Radius',       radius_final_km,           'km';
    'Initial Orbital Velocity',   velocity_initial_orbit,    'm/s';
    'Final Orbital Velocity',     velocity_final_orbit,      'm/s';
    'Semi-major Axis',            semimajor_axis_km,         'km';
    'Periapsis Velocity',         periapsis_velocity,        'm/s';
    'Apoapsis Velocity',          apoapsis_velocity,         'm/s';
    'Initial Orbit Time',         time_initial_hr,           'hour(s)';
    'Final Orbit Time',           time_final_hr,             'hour(s)';
    'Total Transfer Time (s)',    total_orbit_period,        'sec';
    'Total Transfer Time (h)',    total_period_hr,           'hour(s)';
    'Delta V for First Transfer', deltaV_1,                  'm/s';
    'Delta V for Second Transfer',deltaV_2,                  'm/s';
    'Total Delta V',              deltaV_total,              'm/s';
    'Transfer Time (h)',          transfer_time_hr,          'sec';
    'RAAN',                       raan_deg,                  'degrees';
    'Phase Angle',                Phase_AngleD,              'degrees';
};

% Display using a loop
fprintf('\n--------------------------------------------\n');
fprintf('<strong>Hooman Transfer Analysis Results\n</strong>');
fprintf('--------------------------------------------\n');
    for k = 1:size(labels,1)
        fprintf('%s: %.4f %s\n', labels{k,1}, labels{k,2}, labels{k,3});
    end

    if efficiency == 1
        fprintf(2,'\nThe Hohmann method is not efficient for this kind of transfer,\n');
        fprintf(2,'use Bi-Elliptic instead.\n');
    end

    if Phase_AngleD < 0
        fprintf(2,'\nSince the result is negative, you need the Target spacecraft');
        fprintf(2, '\nto be ahead of you by %.4f degrees.\n', Phase_AngleD);
    end

%% FIGURE WITH ANIMATION WILL BE ADDED IN THE LATER VERSION

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