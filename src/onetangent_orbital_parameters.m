function [radius_initial, radius_final, velocity_initial_orbit, velocity_final_orbit, ...
    semimajor_axis, periapsis_velocity, apoapsis_velocity, ...
    time_final_orbit, time_initial_orbit, total_orbit_period, ...
    deltaV_1, deltaV_2, deltaV_total, transfer_time, efficiency, Phase_AngleD] = onetangent_orbital_parameters(initial_altitude, final_altitude, efficiency)
    
    cvc = contantsvalues_convertions(); % Load the constants
    
    % Convert altitudes from kilometers to meters
    radius_initial = (cvc.Earth_radius + initial_altitude) * 1000; % meters
    radius_final   = (cvc.Earth_radius + final_altitude) * 1000;   % meters

    ratio_of_orbital_radius = radius_final / radius_initial; % ratio of final to initial orbital radius
        if ratio_of_orbital_radius < 1
            error('Final altitude must be greater than initial altitude.');
        else if ratio_of_orbital_radius >= 11.94
            efficiency = 1;
        else
            efficiency = 0;
        end

    % Calculate the initial and final velocities based on the gravitational constant
    % the formula for circular orbit velocity: v = sqrt(GM/r)
    % where G is the gravitational constant and M is the mass of the Earth.
    % Here we assume a circular orbit since it is our target orbit for the Hohmann transfer.
    velocity_initial_orbit = sqrt(cvc.Earth_gravitational_constant / radius_initial); % m/s
    velocity_final_orbit   = sqrt(cvc.Earth_gravitational_constant / radius_final);   % m/s

    % Calculate the semimajor axis, the semimajor axis is the average of the initial and final radii
    semimajor_axis = (radius_initial + radius_final) / 2; % meters

    % Calculate the Periapsis and Apoapsis Velocities
    periapsis_velocity = sqrt((cvc.Earth_gravitational_constant) * (2/radius_initial - 1/radius_final)); % m/s
    apoapsis_velocity  = sqrt((cvc.Earth_gravitational_constant) * (2/radius_final - 1/semimajor_axis)); % m/s

    % Ensure the velocities are non-negative
    periapsis_velocity = max(periapsis_velocity, 0);
    apoapsis_velocity  = max(apoapsis_velocity, 0);

    % Compute for individual and total orbital periods
    time_initial_orbit = sqrt((4 * (pi^2) * (radius_initial^3))/cvc.Earth_gravitational_constant);
    time_final_orbit   = sqrt((4 * (pi^2) * (radius_final^3))/cvc.Earth_gravitational_constant);
    total_orbit_period = time_final_orbit + time_initial_orbit;
    transfer_time      = pi * sqrt((semimajor_axis^3)/cvc.Earth_gravitational_constant); % Time for half the transfer orbit

    Phase_Angle  = pi * (1 - sqrt((2 * radius_final) / (radius_initial + radius_final))); % Phase angle in radians
    Phase_AngleD = Phase_Angle * cvc.radians_to_degrees; % Phase angle in degrees

    % Calculate change of velocities of the transfers
    deltaV_1     = periapsis_velocity - velocity_initial_orbit; % m/s
    deltaV_2     = velocity_final_orbit - apoapsis_velocity;   % m/s
    deltaV_total = deltaV_1 + deltaV_2; % m/s
end