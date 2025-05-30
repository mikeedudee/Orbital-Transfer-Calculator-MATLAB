function [initial_orbit, first_half, second_half, final_orbit] = orbits(raan, i_rad, omega, ...
    radius_initial, radius_final, theta, incl, ...
    semimajor_axis_1, eccentricity_1, ...
    semimajor_axis_2, eccentricity_2, use_bielliptic)

    % Initial Circular Orbit
    r0 = radius_initial;
    x0 = r0 * cos(theta);
    y0 = r0 * sin(theta);
    z0 = zeros(size(theta));

    % Final Circular Orbit
    rf = radius_final;
    xf = rf * cos(theta);
    yf = rf * sin(theta);
    zf = zeros(size(theta));

    % Rotation for Inclination and RAAN
    Rz_raan = [cosd(raan) -sind(raan) 0;
               sind(raan)  cosd(raan) 0;
               0           0          1];
    Rx_inc = [1 0 0;
              0 cosd(incl) -sind(incl);
              0 sind(incl)  cosd(incl)];
    R_initial = Rz_raan * Rx_inc;

    initial_orbit = R_initial * [x0; y0; z0];
    final_orbit   = R_initial * [xf; yf; zf];

    if use_bielliptic
        % First Transfer Orbit: Initial to Intermediate
        a1 = semimajor_axis_1;
        e1 = eccentricity_1;
        r1 = a1 * (1 - e1^2) ./ (1 + e1 * cos(theta));
        xt1 = r1 .* cos(theta);
        yt1 = r1 .* sin(theta);
        zt1 = zeros(size(theta));
        transfer_1 = R_initial * [xt1; yt1; zt1];

        % Second Transfer Orbit: Intermediate to Final
        a2 = semimajor_axis_2;
        e2 = eccentricity_2;
        r2 = a2 * (1 - e2^2) ./ (1 + e2 * cos(theta));
        xt2 = r2 .* cos(theta);
        yt2 = r2 .* sin(theta);
        zt2 = zeros(size(theta));
        transfer_2 = R_initial * [xt2; yt2; zt2];

        % Combine both transfer arcs
        half_idx = floor(length(theta)/2);
        first_half  = transfer_1(:, 1:half_idx);
        second_half = transfer_2(:, half_idx:end);
    else
        % Hohmann Transfer Orbit
        a = semimajor_axis_1;
        e = eccentricity_1;
        r = a * (1 - e^2) ./ (1 + e * cos(theta));
        xt = r .* cos(theta);
        yt = r .* sin(theta);
        zt = zeros(size(theta));
        transfer = R_initial * [xt; yt; zt];

        % Divide into two halves
        half_idx = floor(length(theta)/2);
        first_half  = transfer(:, 1:half_idx);
        second_half = transfer(:, half_idx:end);
    end
end
