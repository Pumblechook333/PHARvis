% This script serves to illustrate 2D raytracing results in a planar 
% slice of the IRI Ionosphere. Adapted from PHaRLAP "ray_test1.mat"

%%
% Clear console output and memory

clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")
clf

%% Constants / Settings

saveplots = 0;          % Toggle saving plots in local repository
showplots = 1;          % Toggle showing plots in seperate windows
debug = 1;              % Toggle for light / heavy settings
if debug
    fprintf(['DEBUG MODE: ON' newline newline])
    elevs = 0:2:30;     % Range of elevations for rays
    hours = 15:1:16;      % number of hours (from zero) to generate
    nhops = 1:1:3;      % number of hops to raytrace
else
    fprintf(['DEBUG MODE: OFF' newline newline])
    elevs = 0:0.2:90;
    hours = 0:1:23;       % number of hours (from zero) to generate
    nhops = 1:1:3;                   % number of hops to raytrace
end

UT = [2021 7 1 0 0];                % UT - year, month, day, hour, minute
R12 = 57;                           % R12 index
speed_of_light = 2.99792458e8;      % Constant physical property

honing = 1;                         % Option to only display rays within 
                                    % specified range
RX_range = 2600;                    % Distance of honing target from origin (km)
range = 100;                        % Allowable distance from target (km)
wgs84 = wgs84Ellipsoid('km');       % wgs84 Earth model

num_elevs = length(elevs);          % number of rays 
freq = 10.0;                        % ray frequency (MHz)
freqs = freq.*ones(size(elevs));    % creates array of selected freq

%% Ionosphere settings

TX_coord = [40.67583063, -105.038933178];   % TX Station
origin_lat = TX_coord(1);                   % TX latitude (°)
origin_lon = TX_coord(2);                   % TX longitude (°)
origin_height = 0;                          % Origin height above sea (km)

RX_coord = [40.742018, -74.178975];         % RX Station
receiver_lat = RX_coord(1);                 % latitude of the start point of ray
receiver_lon = RX_coord(2);                 % longitude of the start point of ray
target_height = 250e3;                      % Target height above sea (km)

aim = 0;
ray_bear = bearing(origin_lat, origin_lon, receiver_lat, receiver_lon, ...
                   origin_height, target_height, aim, wgs84);

%% Raytracing Settings

tol = [1e-7 0.01 25];        % ODE tolerance and min/max step sizes
doppler_flag = 1;            % generate ionosphere 5 minutes later so that
                             % Doppler shift can be calculated
irregs_flag = 0;             % no irregularities - not interested in 
                             % Doppler spread or field aligned irregularities
kp = 0;                      % kp not used as irregs_flag = 0. Set it to a 
                             % dummy value 

fprintf( ['\n' ...
  'Example of 2D numerical raytracing for a fan of rays for a' ...
  'WGS84 ellipsoidal Earth\n\n'])

% generate ionospheric, geomagnetic and irregularity grids
max_range = 10000;      % maximum range for sampling the ionosphere (km)
num_range = 201;        % number of ranges (must be < 2000)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)

start_height = 0 ;      % start height for ionospheric grid (km)
height_inc = 3;         % height increment (km)
num_heights = 200;      % number of  heights (must be < 2000)

colors = ["#F05039" "#EEBAB4" "#3D65A5"]; % Colors for different hops

gen = 0;                % Toggle 2D ionosphere generate / load

%%
% Repeat raytracing / optional ionosphere generation per hour and plot

for i = hours
    UT(4) = i;                 % Update the UTC hour
    fprintf('\nHour %d:\n', i);

    % Either Generate or Load 2D ionosphere
    dirname = "IONO/2d_ionosphere/";    % Directory of saved ionospheres
    IONOPATH = dirname + i + ...        % Path to previously saved ionosphere
               "_ionosphere_" + num2str(height_inc * num_heights)...
               + "km_series_" + R12;
    if gen
        tic
        fprintf('Generating ionospheric grid... ')
        [iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
            gen_iono_grid_2d(origin_lat, origin_lon, R12, UT, ray_bear, ...
                             max_range, num_range, range_inc, start_height, ...
                     height_inc, num_heights, kp, doppler_flag, 'iri2020');
        toc
        
        if not(isfolder(dirname))
                    mkdir(dirname)
        end
        
        save(IONOPATH, 'iono_pf_grid', 'iono_pf_grid_5' ,...
            'collision_freq', 'irreg')
    else
        fprintf('Loading ionospheric grid... ')
        load = matfile(IONOPATH);
        iono_pf_grid = load.iono_pf_grid;
        iono_pf_grid_5 = load.iono_pf_grid_5;
        collision_freq = load.collision_freq;
        irreg = load.irreg;
    end

    % convert plasma frequency grid to  electron density in electrons/cm^3
    iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
    iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;
    
    flags = [0 0 0];
    in_range_rays = struct('initial_elev', {}, ... 
                  'frequency', {}, ...
                  'ground_range', {}, ...
                  'height', {}, ...
                  'group_range', {}, ...
                  'phase_path', {}, ...
                  'geometric_distance', {}, ...
                  'electron_density', {}, ...
                  'refractive_index', {}, ...
                  'collision_frequency', {}, ...
                  'absorption', {});
    for h = nhops        
        % call raytrace for a fan of rays
        % first call to raytrace so pass in the ionospheric and geomagnetic grids 
        fprintf('Generating %d 2D NRT rays ...', num_elevs);
        tic
        [ray_data, ray_path_data] = ...
           raytrace_2d(origin_lat, origin_lon, elevs, ray_bear, freqs, h, ...
                       tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
                   collision_freq, start_height, height_inc, range_inc, irreg);
        toc;
        
        if honing          
            for r = ray_path_data
                in_range = chk_dist_2d(range, r, RX_range) == 1;
                if in_range
                    in_range_rays(end+1) = r; % Append the current ray to the new structure array
                end
            end
        else
            in_range_rays = ray_path_data;
        end
        
        flags(h) = length(in_range_rays);
        
    end

    % plot the rays and ionosphere
    figure(i+1)

    if ~showplots
        set(gcf,'visible','off')    % Do not display plot (halts loop)
    end
    

    % Plot window string
    UT_str = [num2str(UT(3)) '/' num2str(UT(2)) '/' num2str(UT(1)) '  ' ...
              num2str(UT(4), '%2.2d') ':' num2str(UT(5), '%2.2d') 'UT'];
    freq_str = [num2str(freq) 'MHz'];
    R12_str = num2str(R12);
    lat_str = num2str(origin_lat);
    lon_str = num2str(origin_lon);
    bearing_str = num2str(ray_bear);
    nhop_str = num2str(h);
    fig_str = [UT_str '   ' freq_str '   R12 = ' R12_str '  TX lat = ' lat_str ...
               ', TX lon = ' lon_str ', bearing = ' bearing_str];
    set(gcf, 'name', fig_str)

    start_range = 0;    % x-lim minimum (TX)
    end_range = 3000;   % x-lim maximum (past RX)
    end_range_idx = fix((end_range-start_range) ./ range_inc) + 1;

    start_ht = start_height;    % Bottom of simulated ionosphere (0km)
    start_ht_idx = 1;           
    end_ht = 400;               % Top of simulated ionosphere (400km)
    end_ht_idx = fix(end_ht ./ height_inc) + 1;

    % Slice total grid to plotted section
    iono_pf_subgrid = iono_pf_grid(start_ht_idx:end_ht_idx, 1:end_range_idx);
    % Plot iono subgrid with rays
    [axis_handle, ray_handle] = ...
    plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, range_inc, ...
        start_ht, end_ht, height_inc, in_range_rays, 'color', [1, 1, 0.99], ...
        'linewidth', 2);

    colormap gray
    
    % Resize figure to screen if showing
    if showplots
        sc_sz = get(0, 'Screensize');   % Get screen dimensions
        sc_sz(1) = 50;                  % Bottom left x pix
        sc_sz(2) = 100;                 % Bottom left y pix
        sc_sz(3) = sc_sz(3) / 1.05;     % Width pix
        sc_sz(4) = sc_sz(4) / 1.8;      % Height pix
        set(gcf, 'Position', sc_sz);    % Set to current figure
    end
    
    % Color the traced ray according to its nhops
    if flags(1) > 0
        set(ray_handle(1:flags(1)), 'color', colors(1), 'linestyle', '-'); 
    end
    if flags(2) > flags(1)
        set(ray_handle(flags(1)+1:flags(2)), 'color', colors(2), 'linestyle', '--');
    end
    if flags(3) > flags(2)
        set(ray_handle(flags(2)+1:flags(3)), 'color', colors(3), 'linestyle', '-.');
    end
    
    % Plot a red line at fractions from center vertical
    div = 250;
    
    step = 100 /div;                    % 100km past 2500km
    marker = 4 + step;                  % 4 ticks right of center
    plot_vertical_line(marker, '2600')  % Line at 2600 (k2mff)
    
    step = 50 /div;                     % 50km past 1250km
    marker = -1 + step;                 % 1 tick left of center
    plot_vertical_line(marker, '1300')  % Line at 1300 (k2mff)

    % Create title at appropriate distance above plot
    plot_title(fig_str)

    % Set the colorbar axes limits
    ax = gca; 
    ax.CLim = [0,8]; 
    ax.Colorbar.Position(2) = 65;

    % Path to 2d raytrace plots
    dirname = "PLOTS/2d_raytrace_combined_cb/";
    filename = dirname + "Z" + UT(4) + "_" + UT(5)+ ".jpg";
    if saveplots
        set(gcf, 'paperorientation', 'portrait')
        % Adjust paperposition to change size of printed plot: 
        % [left, bottom, width, height]
        set(gcf, 'paperunits', 'cent', 'paperposition', [0 0 50 20])
        set(gcf, 'papertype', 'a4') 
        
        if not(isfolder(dirname))
                    mkdir(dirname)
        end
        
        print(filename, '-djpeg', '-r300');   % Save the figure to specified path
    end
end

%% 
% Methods for editing plot_ray_iono_slice output 
% and checking 2D RT reception 

function plot_vertical_line(i, range_str)
% This method will plot a red line at the ticks of the plot.
% 0 is the center tick (not 0km), and the others are accessed by 
% - or + integers. Decimal input allows for a marker to be placed 
% inbetween ticks.

rad_earth = 6371;
start_height = 0;
start_range = 0;
end_range = 3000;
one_tick = 3.14159265 / 16 / 5;
x = one_tick * i;

max_range = end_range - start_range;
tick_len = (max_range ./ 30000) .* 200 + 1000;
tick_r = rad_earth + start_height;

tick_len = 400;
tick_X1 = tick_r .* sin(x);
tick_X2 = (tick_r + tick_len) .* sin(x);
tick_Y1 = tick_r .* cos(x);
tick_Y2 = (tick_r + tick_len) .* cos(x);

adj = 30;

hold on
plot([tick_X1 tick_X2], [tick_Y1 tick_Y2], 'r', 'LineWidth', 3)
text(tick_X2, tick_Y2 + adj, range_str, 'fontsize', 18, ...
     'HorizontalAlignment', 'center', 'Color', 'r')
hold off

end

function plot_title(str)
% This method will plot a red line at the ticks of the plot.
% 0 is the center tick (not 0km), and the others are accessed by 
% - or + integers. Decimal input allows for a marker to be placed 
% inbetween ticks.

rad_earth = 6371;
start_height = 0;
start_range = 0;
end_range = 3000;

max_range = end_range - start_range;
tick_len = (max_range ./ 30000) .* 200 + 1000;
tick_r = rad_earth + start_height;
tick_len = 400;

adj = 130;
text_height = tick_r + tick_len + adj;

title_str = "2D IRI-modeled Ionospheric Raytracing" + newline + str;
hold on
text(0, text_height, title_str, 'fontsize', 18, ...
     'HorizontalAlignment', 'center', ...
     'FontWeight', 'bold')
hold off

end

function ray_bear = bearing(origin_lat, origin_lon, receiver_lat, ...
                            receiver_lon, origin_height, ...
                            target_height, aim, wgs) 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Determine Bearing

    arguments
        origin_lat                  % latitude of TX (°)
        origin_lon                  % longitude of TX (°)
        receiver_lat                % latitude of RX (°)
        receiver_lon                % longitude of RX (°)
        origin_height = 0           % TX height above sea (km)
        target_height = 250e3       % RX height above sea (km)
        aim = 0                     % Adjusted bearing (°)
        wgs = wgs84Ellipsoid        % Spherical earth geometry
    end

    % find the midpoint geographic location between the transmitter and
    % receiver along great circle path
    lat1=origin_lat;
    lon1=origin_lon;
    lat2=receiver_lat;
    lon2=receiver_lon;
    lat1_rad = deg2rad(lat1);
    lon1_rad = deg2rad(lon1);
    lat2_rad = deg2rad(lat2);
    lon2_rad = deg2rad(lon2);

    % find midpoint latitude
    X = cos(lat2_rad) * cos(lon2_rad - lon1_rad);
    Y = cos(lat2_rad) * sin(lon2_rad - lon1_rad);
    mid_lat_rad = atan2(sin(lat1_rad) + sin(lat2_rad), sqrt((cos(lat1_rad) + X) ^ 2 + Y ^ 2));
    % find midpoint longitude
    mid_lon_rad = lon1_rad + atan2(Y, cos(lat1_rad) + X);
    % convert rad into degree
    mid_lat = rad2deg(mid_lat_rad);
    mid_lon = rad2deg(mid_lon_rad);

    [az,~,~] = geodetic2aer(mid_lat,mid_lon,target_height,origin_lat,...
                            origin_lon,origin_height,wgs);
    ray_bear = az + aim;          
    fprintf("Bearing of: " + az + "° \n");

end

function in_range = chk_dist_2d(range, ray_N, RX_range)
                        
            if ~exist('RX_range', 'var')
                RX_range = 2600;
            end
            
            hrange = range;
            
            if ray_N.height(end) < hrange
                ray_range = ray_N.ground_range(end);

                d = abs(RX_range - ray_range);

                % Check if endpoint of ray is within ground range
                if d < range
                    in_range = 1;
                else
                    in_range = 0;
                end

            else
                in_range = -1;
            end
            
end
        