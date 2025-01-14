% This script serves to generate 2D raytracing in a planar slice of the IRI
% Ionosphere.
clc;
clear;

UT = [2021 7 1 0 0];           % UT - year, month, day, hour, minute
R12 = 57;                       % R12 index
speed_of_light = 2.99792458e8;

honing = 1;             % Option to only display rays within 
                        % specified range

debug = 0;              % Toggle for light / heavy settings
if debug
    elevs = 0:2:30;     % Range of elevations for rays
    hours = 0:1:0;       % number of hours (from zero) to generate
    nhops = 1:1:1;                   % number of hops to raytrace
else
    elevs = 0:0.2:90;
    hours = 0:1:23;       % number of hours (from zero) to generate
    nhops = 1:1:3;                   % number of hops to raytrace

end

num_elevs = length(elevs);      % number of rays 
freq = 10.0;                    % ray frequency (MHz)
freqs = freq.*ones(size(elevs));

mode = 1;   % O, No-B, X Mode
gen = 0;    % Toggle generation of 3D Ionosphere
brk = 0;    % Initialize "ray breakdown"
obj = IONS(UT, elevs, freq, R12, mode, gen, brk);   % Load IONS object

ray_bear = obj.ray_bears(1);            % bearing of rays

TX_coord = [40.67583063, -105.038933178];   % WWV Station
origin_lat = TX_coord(1);          % latitude of the start point of ray
origin_long = TX_coord(2);         % longitude of the start point of ray

% tol = [1e-7 .01 10];         % ODE tolerance and min/max step sizes
tol = [1e-7 0.01 25];
doppler_flag = 1;            % generate ionosphere 5 minutes later so that
                             % Doppler shift can be calculated
irregs_flag = 0;             % no irregularities - not interested in 
                             % Doppler spread or field aligned irregularities
kp = 0;                      % kp not used as irregs_flag = 0. Set it to a 
                             % dummy value 

fprintf( ['\n' ...
  'Example of 2D numerical raytracing for a fan of rays for a' ...
  'WGS84 ellipsoidal Earth\n\n'])

%
% generate ionospheric, geomagnetic and irregularity grids
%
max_range = 10000;      % maximum range for sampling the ionosphere (km)
num_range = 201;        % number of ranges (must be < 2000)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)

start_height = 0 ;      % start height for ionospheric grid (km)
height_inc = 3;         % height increment (km)
num_heights = 200;      % number of  heights (must be < 2000)

colors = ['r' 'g' 'b'];
gen = 0;    % Toggle 2D ionosphere generate / load

for h = nhops
    for i = hours
        UT = [2021 7 1 i 0];           % UT - year, month, day, hour, minute
        fprintf('\nHour %d:\n', i);

        % Either Generate or Load 2D ionosphere
        IONOPATH = "IONO/2d_ionosphere/" + i + "_ionosphere_600km_series_" + R12;
        if gen
            tic
            fprintf('Generating ionospheric grid... ')
            [iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
                gen_iono_grid_2d(origin_lat, origin_long, R12, UT, ray_bear, ...
                                 max_range, num_range, range_inc, start_height, ...
                         height_inc, num_heights, kp, doppler_flag, 'iri2020');
            toc

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

        %
        % Example 1 - Fan of rays, 10 MHz, single hop. Print to encapsulated
        % postscript and PNG. Note the transition from E-low to E-High to F2-low modes.
        %

        % call raytrace for a fan of rays
        % first call to raytrace so pass in the ionospheric and geomagnetic grids 
        fprintf('Generating %d 2D NRT rays ...', num_elevs);
        tic
        [ray_data, ray_path_data] = ...
           raytrace_2d(origin_lat, origin_long, elevs, ray_bear, freqs, h, ...
                       tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
                   collision_freq, start_height, height_inc, range_inc, irreg);
        toc;
        
        
        if honing
            RX_range = 2600;
            range = 100;
            wgs84 = wgs84Ellipsoid('km');
            
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
            for r = ray_path_data
                in_range = obj.chk_dist_2d(range, r, RX_range) == 1;
                if in_range
                    in_range_rays(end+1) = r; % Append the current ray to the new structure array
                end
            end
        else
            in_range_rays = ray_path_data
        end

        % plot the rays and ionosphere
        figure(i+1)

        if ~debug
            set(gcf,'visible','off')    % Do not display plot (halts loop)
        end

        % Plot window string
        UT_str = [num2str(UT(3)) '/' num2str(UT(2)) '/' num2str(UT(1)) '  ' ...
                  num2str(UT(4), '%2.2d') ':' num2str(UT(5), '%2.2d') 'UT'];
        freq_str = [num2str(freq) 'MHz'];
        R12_str = num2str(R12);
        lat_str = num2str(origin_lat);
        lon_str = num2str(origin_long);
        bearing_str = num2str(ray_bear);
        nhop_str = num2str(h);
        fig_str = [UT_str '   ' freq_str '   R12 = ' R12_str '  TX lat = ' lat_str ...
                   ', TX lon = ' lon_str ', bearing = ' bearing_str ...
                   ', nhops = ' nhop_str];
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
        plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, range_inc, ...
            start_ht, end_ht, height_inc, in_range_rays, 'color', [1, 1, 0.99], ...
            'linewidth', 2);

        % Plot a red line at fractions from center vertical
        % Line at 2600 (k2mff)
        step = 100/250; % 100km past 2500km
        marker = 4 + step;
        plot_vertical_line(marker, '2600')
        % Line at 1300 (k2mff)
        step = 50/250;  % 50km past 1250km
        marker = -1 + step;
        plot_vertical_line(marker, '1300')

        % Text for plot title
        plot_title(fig_str)

        % Set the colorbar axes limits
        ax = gca; 
        ax.CLim = [0,8]; 
        ax.Colorbar.Position(2) = 65;

    %     set(gcf,'units','normal')
    %     pos = get(gcf,'position');
    % %     pos(2) = 0.55;
    %     pos(2) = 0.65;
    %     set(gcf,'position', pos)

        % uncomment the following to print figure to hi-res ecapsulated postscript
        % and PNG files
        set(gcf, 'paperorientation', 'portrait')
    %     set(gcf, 'paperunits', 'cent', 'paperposition', [0 0 61 18])
        % Adjust paperposition to change size of printed plot: 
        % [left, bottom, width, height]
        set(gcf, 'paperunits', 'cent', 'paperposition', [0 0 50 20])
        set(gcf, 'papertype', 'a4') 
        % print -depsc2 -loose -opengl test.ps 

        % Path to 2d raytrace plots
        filename = "PLOTS/2d_raytrace_honed/" + h + "_hop_" + UT(4) + "_" + UT(5)+ ".png";
        % print -dpng filename
        if ~debug
            print(filename, '-dpng');   % Save the figure to specified path
        end
    end
end

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

