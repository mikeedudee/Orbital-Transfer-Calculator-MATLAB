function figure_main()

    % Figure setup
    x      = 150; % x-coordinate of the figure's bottom-left corner
    y      = 0;   % y-coordinate of the figure's bottom-left corner

    width  = 1080; % width in pixels
    height = 720;  % height in pixels

    fig = figure('Position', [x, y, width, height], 'Toolbar', 'none', 'NumberTitle', 'off', 'Name', 'Hohmann Transfer');
    hold on;
    grid on;

    % Improve realism
    axis equal;
    axis off;
    shading interp;  % Smooth shading
    view(50, 30);    % Adjust camera
    rotate3d on;
    


    %% Custom menu
    % Add read lincense menu item
    license = uimenu('Text', 'Read License');
    set(0, 'showhiddenhandles', 'off');
    LicenseHandles = findall(fig);
    set(LicenseHandles(2), 'foregroundColor', [1 0 0])
    uimenu(license, 'Label', 'Program License', 'Callback', @license_callback); 

    % Add a readme menu item
    readme = uimenu('Text', 'ReadMe');
    set(0, 'showhiddenhandles', 'off');
    ReadMeHandles = findall(fig);
    set(ReadMeHandles(2), 'foregroundColor', [1 0 0])
    uimenu(readme, 'Label', 'ReadMe', 'Callback', @readme_callback); 

    % Add a paper menu item
    paper = uimenu('Text', 'Read the Methodological Paper');
    set(0, 'showhiddenhandles', 'off');
    paperHandles = findall(fig);
    set(ReadMeHandles(2), 'foregroundColor', [0 1 0])
    uimenu(paper, 'Label', 'Link Here', 'Callback', @paper_callback); 

    % Custome menu callback functions
    function hyperlink_callback(source, event)
        web('https://github.com/mikeedudee/Orbital-Transfer-Calculator-MATLAB/blob/32a5e8e1671f3f1d20fbab97c10c6a0e0a213b52/LICENSE', '-browser');
    end

    function readme_callback(source, event)
        web('https://github.com/mikeedudee/Orbital-Transfer-Calculator-MATLAB/blob/32a5e8e1671f3f1d20fbab97c10c6a0e0a213b52/README.md', '-browser');
    end

    function paper_callback(source, event)
        %web('');
    end

end