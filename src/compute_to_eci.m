function compute_to_eci()

    % Loop to compute ECI position for each true anomaly
    for idx = 1:length(theta)
        nu = theta(idx); % true anomaly
        
        % Radius from orbit equation
        r = (semimajor_axis * (1 - eccentricity^2)) / (1 + eccentricity * cos(nu));
        
        % Position vector in PQW frame
        r_pqw = [r * cos(nu); r * sin(nu); 0];
        
        % Rotation matrices for transformation PQW -> ECI
        R3_W = [cos(-RAAN_rad) -sin(-RAAN_rad) 0;
                sin(-RAAN_rad)  cos(-RAAN_rad) 0;
                0               0              1];
        
        R1_i = [1 0 0;
                0 cos(-i_rad) -sin(-i_rad);
                0 sin(-i_rad)  cos(-i_rad)];
            
        R3_w = [cos(-omega) -sin(-omega) 0;
                sin(-omega)  cos(-omega) 0;
                0            0           1];
        
        Q_pqw2eci = R3_W * R1_i * R3_w;
        r_eci(:, idx) = Q_pqw2eci * r_pqw;
    end

end