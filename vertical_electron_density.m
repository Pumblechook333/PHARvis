%% Setup
clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")

clf
%% Constants / Settings
mode_keys = ["O", "X"];
R12_sel = [-1, 25, 50, 100, 200];
date = [2021 7 1 0 0];
el_start = 0;

hi_res = 0;
if hi_res
    el_inc = 0.2;
    el_stop = 50;
else
    el_inc = 5;
    el_stop = 90;
end

elevs = el_start : el_inc : el_stop;

freq = 10;

gen = 0; % 0 = no gen, 1 = gen
brk = false;

%% GET necessary vars

R12 = R12_sel(1);

elevs_string = " || Initial Elevations: " ...
               + el_start + ":" + el_inc + ":" + el_stop;
r12_string = " || R12: " + R12;

mode = 1;
obj = IONS(date, elevs, freq, R12, mode, gen, brk);

% Raytrace rays at all nhops for certain hour
% hour = 1; % 0th hour = 1
% rays_N = obj.get_hour_rays(hour);

nhops = 1;
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes

hour_field = "i1";
iono_pf_grid = obj.iono_series.(hour_field).iono_pf_grid;
iono_pf_grid_5 = obj.iono_series.(hour_field).iono_pf_grid_5;
collision_freq = obj.iono_series.(hour_field).collision_freq;
Bx = obj.iono_series.(hour_field).Bx;
By = obj.iono_series.(hour_field).By;
Bz = obj.iono_series.(hour_field).Bz;

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;
                
elev = 90; % Shoot a ray straight up (Ionosonde)
bear = 0;
freq = 10;

start_height = 0;
m_latitude = 41.75;
m_longitude = -89.62;

[ray_data_N, ray_N, ~] = ...
            raytrace_3d(m_latitude, m_longitude, start_height, elev, bear, freq, ...
            obj.mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
            collision_freq, obj.iono_grid_parms, Bx, By, Bz, ...
            obj.geomag_grid_parms); 
        
can_plot = false;
if ray_data_N.ray_label == -2
    fprintf("Ray Penetrated Ionosphere.\n")
    can_plot = true;
else
    fprintf("Ray HAS NOT Penetrated Ionosphere.\n")
end

if can_plot
    density = ray_N.electron_density;
    altitude = ray_N.height;
    semilogx(density, altitude, "LineWidth", 2);

    height_sel = 200;
    yline(height_sel, "LineWidth", 2, "Color", 'r');

    legend_cells = {"Density Along Vertical",...
                    "Altitude (" + height_sel + ")"};
    % legend(legend_cells, 'Location', 'northwest');
    xlabel('Density (electrons / cm^3)');
    ylabel('Elevation (km)');
    ylim([0,400])
    grid on;

    set(gca,"FontSize",20)

    ti = "IRI Electron Density at [" + m_latitude + "," + (m_longitude+360) + "]";
    title(ti+r12_string)
end
    
fprintf("Done.\n")
