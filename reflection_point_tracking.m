 clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n")

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% GET necessary vars

date = [2021 7 1 0 0];

el_start = 0;
el_inc = 0.2;
% el_inc = 2;
el_stop = 90;
elevs = el_start : el_inc : el_stop;

freq = 10;
R12_sel = [-1, 25, 50, 100, 200];
R12 = R12_sel(1);
gen = 0; % 0 = no gen, 1 = gen

elevs_string = " || Initial Elevations: " ...
               + el_start + ":" + el_inc + ":" + el_stop;
r12_string = " || R12: " + R12;

obj = IONS(date, elevs, freq, R12, gen);

iono_series = obj.get_iono_series();

iono = obj.get_iono_parms();
iono_grid_parms = iono.iono_grid_parms;
geomag_grid_parms = iono.geomag_grid_parms;
iono_height = iono.iono_height;

g = obj.get_gen_params();
elevs = g.elevs;
freqs = g.freqs;
UT = g.UT;

c = obj.get_coords();
origin_lat = c.origin_lat;
origin_lon = c.origin_lon;
origin_ht = c.origin_ht;
RX_coord = c.RX_coord;
TX_coord = c.TX_coord;

ray_bears = obj.get_bearing();

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plotting

rsta = 1;
rinc = 1;
rsto = 24;

tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
num_elevs = length(elevs); 

% OX_mode - polarization mode of rays: 1 = O, -1 = X, 0 = no field
OX_mode = 1;

max_heights = struct();

range = 100;
wgs84 = wgs84Ellipsoid('km');

nhop_max = 4;
for nhops = 1:1:nhop_max
    hop_field = "hop_" + nhops;
    
    height_arr = zeros(rsto,num_elevs);
    for i = rsta:rinc:rsto
        hour_field = 'i' + string(i);

        iono_pf_grid = iono_series.(hour_field).iono_pf_grid;
        iono_pf_grid_5 = iono_series.(hour_field).iono_pf_grid_5;
        collision_freq = iono_series.(hour_field).collision_freq;
        Bx = iono_series.(hour_field).Bx;
        By = iono_series.(hour_field).By;
        Bz = iono_series.(hour_field).Bz;

        % convert plasma frequency grid to  electron density in electrons/cm^3
        iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
        iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

        [ray_data_N, ray_N, ray_sv_N] = ...
          raytrace_3d(origin_lat, origin_lon, origin_ht, elevs, ray_bears, freqs, ...
                      OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
                      collision_freq, iono_grid_parms, Bx, By, Bz, ...
                      geomag_grid_parms);
                          
        heights = zeros(1,num_elevs);
        % Find maximum height of rays
        
        for elev = 1:1:num_elevs
            grounded_ray = ray_data_N(elev).ray_label == 1;
            if grounded_ray
                in_range = IONS.chk_dist(range, ray_N(elev), wgs84) == 1;
                if in_range
                    h = max(ray_N(elev).height);
                    heights(elev) = h;
                end
            end
        end

        height_arr(i,:) = heights;

    end
    max_heights.(hop_field) = height_arr;
    % finished ray tracing with this ionosphere so clear it out of memory
    clear raytrace_3d
end

clf
figure(1)
pos = get(gcf, 'position');
pos(3) = pos(3)*1.5;
pos(4) = pos(4)*1.5;
set(gcf, 'position', pos)
set(gca, 'Ylim', [0 600])

hold on

colors = ['r', 'g', 'b', 'k'];

plot(0,0, 'color', 'r', 'marker', 'x', 'markersize', 1)
plot(0,0, 'color', 'g', 'marker', 'x', 'markersize', 1)
plot(0,0, 'color', 'b', 'marker', 'x', 'markersize', 1)
plot(0,0, 'color', 'k', 'marker', 'x', 'markersize', 1)
legend('1-hop', '2-hop', '3-hop', '4-hop', 'AutoUpdate', 'off')

for nhops = 1:1:nhop_max
    hop_field = "hop_" + nhops;
    fprintf("Plotting " + nhops + "-hop rays \n")
    
    hr_range = 0:1:24;

    for k = rsta:rinc:rsto
        hour_field = 'i' + string(k);
        
        for i = 1:1:num_elevs
            y = max_heights.(hop_field)(k,i);
            in_range = (0 < y) & (y < iono_height);

            if in_range
                plot(k-1, y, 'color', colors(nhops), 'marker', 'x', ...
                    'markersize', 20, 'LineWidth', 4)
            end
        end
        
    end
    
end

hold off
grid on
fs = 30;

ti = "Maximum Height of Rays across Day";
title((ti + r12_string + elevs_string), 'FontSize', fs)
xlabel('Time (UT)', 'FontSize', fs)
ylabel('Height (km)', 'FontSize', fs)

xticks(0:1:24)

ax = gca;
ax.FontSize = fs/1.5; 

