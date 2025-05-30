function earth_sphere()

    N                 = 1000; % Resolution of the mesh and graphical representation
    earth_img         = imread('earth.jpg'); % Load Earth image
    earth_img_resized = flipud(imresize(earth_img, [N+1, N+1]));  % Match mesh

    % Generate sphere coordinates
    [sx, sy, sz] = sphere(N);
    earth_sphere = surf(sx, sy, sz);
    set(earth_sphere, ...
        'FaceColor', 'texturemap', ...
        'CData', earth_img_resized, ...
        'EdgeColor', 'none');

end