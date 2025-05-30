function animation(first_half, second_half, final_orbit, initial_orbit, Earth_radius_m, rocket_marker, use_bielliptic)

    %% ANIMATION SECTION
    animate_trajectory_segment(initial_orbit, rocket_marker, Earth_radius_m);  % Initial orbit
    animate_trajectory_segment(first_half, rocket_marker, Earth_radius_m);     % First transfer arc

    if use_bielliptic
        % Animate second elliptical arc for bi-elliptic transfer
        animate_trajectory_segment(second_half, rocket_marker, Earth_radius_m); 
    end

    % Match closest point in final orbit to last transfer segment
    % Determine last transfer point depending on the type of transfer
    if use_bielliptic
        last_point = second_half(:, end);
    else
        last_point = first_half(:, end);
    end

    % Match closest point in final orbit to last transfer segment
    distances  = vecnorm(final_orbit - last_point, 2, 1);
    [~, start_idx] = min(distances);


    animate_trajectory_segment(final_orbit(:, start_idx:end), rocket_marker, Earth_radius_m);

    % Continuous animation of final orbit
    while true
        animate_trajectory_segment(final_orbit, rocket_marker, Earth_radius_m);
    end

    %% Internal animation function
    function animate_trajectory_segment(trajectory, marker_handle, Earth_radius_m)
        x = trajectory(1, :) ./ Earth_radius_m;
        y = trajectory(2, :) ./ Earth_radius_m;
        z = trajectory(3, :) ./ Earth_radius_m;

        for i = 1:length(x)
            set(marker_handle, 'XData', x(i), 'YData', y(i), 'ZData', z(i));
            drawnow;
        end
    end

end
