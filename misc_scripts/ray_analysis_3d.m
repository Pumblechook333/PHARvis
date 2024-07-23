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

nhops = 4;                  % number of hops
elevs = [20:1:30];               % initial elevation of rays
freqs = ones(size(elevs))*10;   % frequency (MHz)
speed_of_light = 2.99792458e8;
R12 = -1;
doppler_flag = 0;               % interested in Doppler shift

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
% Determine Bearing

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

wgs84 = wgs84Ellipsoid;

H = 250e3; % (F-layer)
[az,~,~] = geodetic2aer(mid_lat,mid_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
% aim = -2;   % Adjust bearing
aim = 0;
az = az + aim;
% az = 0;
ray_bears = ones(size(elevs))*az ; % initial bearing of rays

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Ionospheric Grid Parameters

% Bottom left corner of grid (lat, lon) in [deg]
grid_corner = [40, -(105+(48/60))];
% grid_corner = [38, -106];
% grid_corner = [36, -110];

ht_start = 0;          % start height for ionospheric grid (km)
ht_inc = 4;             % height increment (km)
num_ht = 100 +1;     

lat_start = grid_corner(1);
lat_inc = 0.5;
num_lat = 4 +1;

lon_start= grid_corner(2);
lon_inc = 1.0;
num_lon = 40 +1;

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

fname_append = "_400km.mat";
% fname = 'ionosphere_800km.mat';
fname = 'ionosphere' + fname_append;

gen = 0;
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

plot3(ray_O(1).lat, ray_O(1).lon, ray_O(1).height, '.b', ...
      'markersize', 5)
set(gca, 'Zlim', [0 400])
hold on
plot3(ray_X(1).lat,  ray_X(1).lon, ray_X(1).height, '.r',  ...
      'markersize',5)
plot3(ray_N(1).lat,  ray_N(1).lon, ray_N(1).height, 'g')
for ii = 2:1:num_elevs
  plot3(ray_O(ii).lat, ray_O(ii).lon, ray_O(ii).height, '.b', ...
        'markersize', 5)
  plot3(ray_X(ii).lat, ray_X(ii).lon, ray_X(ii).height, '.r', ...
        'markersize', 5)
  plot3(ray_N(ii).lat,  ray_N(ii).lon, ray_N(ii).height, 'g')
end  

TX = round([TX_coord(1), TX_coord(2), 0], 2);

% BP = [mid_lat, mid_lon+360, 200];
BP = round([42.12, -82.60, 200], 2);

RX = round([40.79, -69.50, 0], 2);

points = [TX; BP; RX];
fs = 30;

for i = 1:1:3
    x = points(i,1); y = points(i,2); z = points(i,3);
    txt = "[" + x + ", " + y + ", " + z + "]";
    plot3(x, y, z, 'x', 'color', 'k', 'markersize', 20, 'LineWidth',4)
    text(x - 0.5, y, z, txt, 'FontSize', fs-10, 'FontWeight', 'bold') 
end

RX = round([RX_coord(1), RX_coord(2)+360, 0], 2);
x = RX_coord(1); y = RX_coord(2); z = 0;
txt = "[" + x + ", " + y + ", " + z + "]";
plot3(x, y, z, 'o', 'color', 'y', 'markersize', 20, 'LineWidth',4)
text(x - 1, y - 4, z, txt, 'FontSize', fs-10, 'FontWeight', 'bold') 

set(gca, 'XDir','reverse')
hold off
grid on
ti = "3D Raytracing Rays from WWV to K2MFF";
su = UT(1) + "-" + UT(2) + "-" + UT(3) + " | " + ...
           "TX: [" + TX(1) + ", " + TX(2) + "] --> RX: [" +...
           round(RX_coord(1), 2) + ", " + round(RX_coord(2), 2) + "] | " + ...
           "Midpoint : [" + BP(1) + ", " + BP(2) + "]";
       
% title([ti,su], 'FontSize', fs)
xlabel('Latitude (°)', 'FontSize', fs)
ylabel('Longitude (°)', 'FontSize', fs)
zlabel('Height (km)', 'FontSize', fs)
legend('O Mode', 'X Mode', 'No Mag-field', 'FontSize', fs)

ax = gca;
ax.FontSize = fs/1.5; 
