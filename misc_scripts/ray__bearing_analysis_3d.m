clear
clc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% General Parameters

tsel = 1;
switch tsel
    case 1
        UT = [2021 7 1 1 34];           % Sunset
        suntime = "sunset";
    case 2
        UT = [2021 7 1 10 32];        % Sunrise
        suntime = "sunrise";
    case 3
        UT = [2021 7 1 18 03];        % Solar Noon
        suntime = "solar noon";
end

nhops = 10;                         % number of hops
ray_bears = [0:2:360] ;             % initial bearing of rays
elevs = ones(size(ray_bears))*1;    % initial elevation of rays
freqs = ones(size(ray_bears))*10;   % frequency (MHz)
speed_of_light = 2.99792458e8;
R12 = -1;
doppler_flag = 1;               % interested in Doppler shift

fprintf(nhops + "-hop transmission at " + suntime + " \n");

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coordinates of Towers

TX_coord = [40.67583063, -105.038933178];   % WWV Station
RX_coord = [40.742018, -74.178975];         % K2MFF Station
origin_lat = TX_coord(1);             % latitude of the start point of rays
origin_lon = TX_coord(2);            % longitude of the start point of rays
origin_ht = 0.0;                % altitude of the start point of rays
receiver_lat = RX_coord(1); 
receiver_lon = RX_coord(2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Ionospheric Grid Parameters

% Bottom left corner of grid (lat, lon) in [deg]
% grid_corner = [40, -(105+(48/60))];
% grid_corner = [38, -106];
% grid_corner = [36, -110];
grid_corner = [31, -118];


ht_start = 0;          % start height for ionospheric grid (km)
ht_inc = 4;             % height increment (km)
num_ht = 100 +1;     

lat_start = grid_corner(1);
lat_inc = 1;
num_lat = 16 +1;

lon_start= grid_corner(2);
lon_inc = 1.0;
num_lon = 60 +1;

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Magnetospheric Grid Parameters

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = ht_inc;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);

B_lat_start = lat_start;
B_lat_inc = lat_inc;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);

B_lon_start = lon_start;
B_lon_inc = lon_inc;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc); 

geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate Ionosphere

% fname_append = "_400km.mat";
fname_append = "_400km_US.mat";
% fname = 'ionosphere_800km.mat';
fname = 'ionosphere' + fname_append;

gen = 1;
if gen == 1
    fprintf('Generating ionospheric and geomag grids... \n')
    
    tic
    
    [iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, doppler_flag);
    
    save(fname, 'iono_pf_grid', 'iono_pf_grid_5' ,...
        'collision_freq', 'Bx', 'By', 'Bz')
    
    toc
else
    fprintf('Loading ionospheric and geomag grids... \n')
    load = matfile(fname);
    iono_pf_grid = load.iono_pf_grid;
    iono_pf_grid_5 = load.iono_pf_grid_5;
    collision_freq = load.collision_freq;
    Bx = load.Bx;
    By = load.By;
    Bz = load.Bz;
end

fprintf('\n')

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

%
% call raytrace
%
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
num_elevs = length(elevs);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plotting

% Generate the O mode rays
OX_mode = 1;
 
[ray_data_O, ray_O, ray_state_vec_O] = ...
  raytrace_3d(origin_lat, origin_lon, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms);

for rayId=1:num_elevs
  num = length(ray_O(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_O(rayId).lat;
  lon = ray_O(rayId).lon; 
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_lon,'wgs84')/1000.0;
  ray_O(rayId).ground_range = ground_range;
end


% Generate the X mode rays - note in the raytrace_3d call the ionosphere does
% not need to be passed in again as it is already in memory
OX_mode = -1;
 
[ray_data_X, ray_X, ray_sv_X] = ...
  raytrace_3d(origin_lat, origin_lon, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol);

for rayId=1:num_elevs
  num = length(ray_X(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_X(rayId).lat;
  lon = ray_X(rayId).lon;    
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_lon,'wgs84')/1000.0;
  ray_X(rayId).ground_range = ground_range;
end


% Generate the rays for the case where the magnetic field is ignored  - note
% in the raytrace_3d call the ionosphere does not need to be passed in again
% as it is already in memory
OX_mode = 0;

[ray_data_N, ray_N, ray_sv_N] = ...
  raytrace_3d(origin_lat, origin_lon, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol);

for rayId=1:num_elevs
  num = length(ray_N(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_N(rayId).lat;
  lon = ray_N(rayId).lon;    
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_lon,'wgs84')/1000.0;
  ray_N(rayId).ground_range = ground_range;
end

fprintf('\n')


% finished ray tracing with this ionosphere so clear it out of memory
clear raytrace_3d


% plot the rays
figure(1)
pos = get(gcf, 'position');
pos(3) = pos(3)*1.5;
pos(4) = pos(4)*1.5;
set(gcf, 'position', pos)
plot3(ray_O(1).lat, mod(ray_O(1).lon, 360), ray_O(1).height, '.b', ...
      'markersize', 5)
set(gca, 'Zlim', [0 500])
hold on
plot3(ray_X(1).lat,  mod(ray_X(1).lon, 360), ray_X(1).height, '.r',  ...
      'markersize',5)
plot3(ray_N(1).lat,  mod(ray_N(1).lon, 360), ray_N(1).height, 'g')
for ii = 3:2:num_elevs
  plot3(ray_O(ii).lat, mod(ray_O(ii).lon, 360), ray_O(ii).height, '.b', ...
        'markersize', 5)
  plot3(ray_X(ii).lat, mod(ray_X(ii).lon, 360), ray_X(ii).height, '.r', ...
        'markersize', 5)
  plot3(ray_N(ii).lat,  mod(ray_N(ii).lon, 360), ray_N(ii).height, 'g')
end  
set(gca, 'XDir','reverse')
hold off
grid on
xlabel('latitude (deg)')
ylabel('longitude (deg)')
zlabel('Height (km)')
legend('O Mode', 'X Mode', 'No Mag-field')

figure(2)
start_range = 0;
end_range = 3200;
range_inc = 50;
end_range_idx = fix((end_range-start_range) ./ range_inc) + 1;
start_ht = 0;
start_ht_idx = 1;
height_inc = 5;
end_ht = 350;
end_ht_idx = fix(end_ht ./ height_inc) + 1;
iono_pf_subgrid = zeros(end_ht_idx, end_range_idx);
plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, range_inc, ...
    start_ht, end_ht, height_inc, ray_O, 'color', [1, 1, 0.99], 'linewidth', 1);
set(findall(gcf,'-property','FontSize'),'FontSize',10)

