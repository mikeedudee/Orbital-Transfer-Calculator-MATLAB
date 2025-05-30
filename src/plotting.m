function [rocket_marker] = plotting(Earth_radius_m, initial_orbit, first_half, second_half, final_orbit, is_bielliptic)

    % Initialize rocket marker at the first point of the initial orbit
    rocket_marker = plot3(initial_orbit(1,1) ./ Earth_radius_m, ...
                          initial_orbit(2,1) ./ Earth_radius_m, ...
                          initial_orbit(3,1) ./ Earth_radius_m, ...
                          'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 5);

    % Plot initial circular orbit
    plot3(initial_orbit(1,:) ./ Earth_radius_m, ...
          initial_orbit(2,:) ./ Earth_radius_m, ...
          initial_orbit(3,:) ./ Earth_radius_m, 'g', 'LineWidth', 1);

    % Plot first half of the transfer orbit (always solid)
    plot3(first_half(1,:) ./ Earth_radius_m, ...
          first_half(2,:) ./ Earth_radius_m, ...
          first_half(3,:) ./ Earth_radius_m, 'r-', 'LineWidth', 1);

    % Plot second half of the transfer orbit
    if is_bielliptic
        line_style = 'r-';  % Solid line for bi-elliptic
    else
        line_style = 'r--'; % Dashed line for Hohmann
    end

    plot3(second_half(1,:) ./ Earth_radius_m, ...
          second_half(2,:) ./ Earth_radius_m, ...
          second_half(3,:) ./ Earth_radius_m, line_style, 'LineWidth', 1);

    % Plot final orbit
    plot3(final_orbit(1,:) ./ Earth_radius_m, ...
          final_orbit(2,:) ./ Earth_radius_m, ...
          final_orbit(3,:) ./ Earth_radius_m, 'c', 'LineWidth', 1);

    % Small sphere for rocket visualization
    xr = initial_orbit(1,1) / Earth_radius_m;
    yr = initial_orbit(2,1) / Earth_radius_m;
    zr = initial_orbit(3,1) / Earth_radius_m;

    [rx, ry, rz] = sphere(8);
    rocket_size = 0.02;
    surf(rocket_size*rx + xr, rocket_size*ry + yr, rocket_size*rz + zr, ...
         'FaceColor', 'r', 'EdgeColor', 'none');

end
