function [radius_initial, radius_final, velocity_initial_orbit, velocity_final_orbit, ...
    semi_major_axis_first, eccentricity_first, periapsis_velocity_i, apoapsis_velocity_i, ...
    semi_major_axis_second, eccentricity_second, periapsis_velocity_ii, apoapsis_velocity_ii, ...
    time_initial_orbit, time_intermediate_orbit, time_final_orbit, ...
    total_orbit_period, transfer_time_first, transfer_time_second, total_transfer_time, ...
    deltaV_1, deltaV_2, deltaV_3, deltaV_total, efficiency, radius_intermediate] = hohmann_orbital_parameters(initial_altitude, final_altitude, alt_intermediate, inclination_c, efficiency)

    cvc = contantsvalues_convertions(); % Load constants

    % Convert altitudes from kilometers to meters
    radius_initial      = (cvc.Earth_radius + initial_altitude) * 1000; % meters
    radius_final        = (cvc.Earth_radius + final_altitude) * 1000;   % meters
    radius_intermediate = (cvc.Earth_radius + alt_intermediate) * 1000;

    % Orbital radius ratio check
    ratio_of_orbital_radius = radius_final / radius_initial;
    if ratio_of_orbital_radius < 1
        error('Final altitude must be greater than initial altitude.');
    elseif ratio_of_orbital_radius <= 11.94
        efficiency = 0;
    else
        efficiency = 1;
    end

    %% First Circular Orbit
    velocity_initial_orbit = sqrt(cvc.Earth_gravitational_constant / radius_initial);

    %% First Transfer Ellipse
    semi_major_axis_first = (radius_initial + radius_intermediate) / 2;
    eccentricity_first    = abs((radius_intermediate - radius_initial) / ...
                              (radius_intermediate + radius_initial)); % non-negative

    periapsis_velocity_i = velocity_initial_orbit * sqrt(1 + eccentricity_first * cos(inclination_c));
    apoapsis_velocity_i  = sqrt(cvc.Earth_gravitational_constant / radius_intermediate) * ...
                           sqrt(1 - eccentricity_first * cos(inclination_c));

    %% Second Transfer Ellipse
    semi_major_axis_second = (radius_intermediate + radius_final) / 2;
    eccentricity_second    = abs((radius_intermediate - radius_final) / ...
                              (radius_final + radius_intermediate)); % non-negative

    periapsis_velocity_ii = sqrt(cvc.Earth_gravitational_constant / radius_final) * ...
                            sqrt(1 + eccentricity_second * cos(inclination_c));
    apoapsis_velocity_ii  = sqrt(cvc.Earth_gravitational_constant / radius_intermediate) * ...
                            sqrt(1 - eccentricity_second * cos(inclination_c));

    % Ensure velocities are non-negative
    periapsis_velocity_ii = max(periapsis_velocity_ii, 0);
    apoapsis_velocity_ii  = max(apoapsis_velocity_ii, 0);

    %% Final Circular Orbit
    velocity_final_orbit = sqrt(cvc.Earth_gravitational_constant / radius_final);

    %% Orbital Times
    time_initial_orbit      = 2 * pi * sqrt(radius_initial^3 / cvc.Earth_gravitational_constant);
    time_intermediate_orbit = 2 * pi * sqrt(radius_intermediate^3 / cvc.Earth_gravitational_constant);
    time_final_orbit        = 2 * pi * sqrt(radius_final^3 / cvc.Earth_gravitational_constant);

    total_orbit_period      = time_initial_orbit + time_intermediate_orbit + time_final_orbit;
    transfer_time_first     = pi * sqrt(semi_major_axis_first^3 / cvc.Earth_gravitational_constant);
    transfer_time_second    = pi * sqrt(semi_major_axis_second^3 / cvc.Earth_gravitational_constant);
    total_transfer_time     = transfer_time_first + transfer_time_second;

    %% Delta V Calculations
    deltaV_1 = periapsis_velocity_i - velocity_initial_orbit; % m/s
    deltaV_2 = apoapsis_velocity_ii - apoapsis_velocity_i;    % m/s
    deltaV_3 = periapsis_velocity_ii - velocity_final_orbit; % m/s
    deltaV_total = deltaV_1 + deltaV_2 + deltaV_3; % m/s 
    
end
