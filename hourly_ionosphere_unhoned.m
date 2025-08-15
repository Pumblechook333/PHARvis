% This is the script which produced the 3D figure that is
% currently visible in the paper, as of 1/14/25 at
% 2:18PM

clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n")

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% GET necessary vars

% DEBUG
debug = 0;

date = [2021 7 1 0 0];

if debug
    elevs = 0:2:90;
else
    elevs = 0:0.2:90;
end

freq = 10;
R12 = 57;
mode = 1;
gen = 0;
brk = 0;
obj = IONS(date, elevs, freq, R12, mode, gen, brk);

iono_series = obj.get_iono_series();

iono = obj.get_iono_parms();
iono_grid_parms = iono.iono_grid_parms;
geomag_grid_parms = iono.geomag_grid_parms;
iono_height = iono.iono_height;

g = obj.get_gen_params();
elevs = g.elevs;
freqs = g.freqs;
UT = repmat(date, [25,1]);
for i=1:1:25
    UT(i,4) = i-1;
end


c = obj.get_coords();
origin_lat = c.origin_lat;
origin_lon = c.origin_lon;
origin_ht = c.origin_ht;
RX_coord = c.RX_coord;
TX_coord = c.TX_coord;

ray_bears = obj.get_bearing();

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plotting

rsta = 17;
rsto = 20;
nhops = 2;

tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
num_elevs = length(elevs); 
modecol = ['r', 'g', 'b'];
% Ellipsoid to base distance calcualtions on
wgs84 = wgs84Ellipsoid('km');
% Range for endpoint of ray (km)
range = 100;

rot=90;
i = rsta;
while i <= rsto
    field = 'i' + string(i);
    
    iono_pf_grid = iono_series.(field).iono_pf_grid;
    iono_pf_grid_5 = iono_series.(field).iono_pf_grid_5;
    collision_freq = iono_series.(field).collision_freq;
    Bx = iono_series.(field).Bx;
    By = iono_series.(field).By;
    Bz = iono_series.(field).Bz;
    
    % convert plasma frequency grid to  electron density in electrons/cm^3
    iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
    iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;
    
    clf
    figure(1)
    pos = get(gcf, 'position');
    pos(3) = pos(3)*1.5;
    pos(4) = pos(4)*1.5;
    set(gcf, 'position', pos)
    set(gca, 'Zlim', [0 iono_height])
    
    plot3(TX_coord(1), TX_coord(2), 0, 'r', 'DisplayName', 'X-mode')
    hold on
    % plot3(TX_coord(1), TX_coord(2), 0, 'g', 'DisplayName', 'No B')
    plot3(TX_coord(1), TX_coord(2), 0, 'b', 'DisplayName', 'O-mode')
    legend('Location', 'best');
    
    for oxmode = -1:2:1
        plotcol = modecol(oxmode + 2);
        
        [ray_data_N, ray_N, ray_sv_N] = ...
          raytrace_3d(origin_lat, origin_lon, origin_ht, elevs, ray_bears, freqs, ...
                      oxmode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
                      collision_freq, iono_grid_parms, Bx, By, Bz, ...
                      geomag_grid_parms);

        % finished ray tracing with this ionosphere so clear it out of memory
        clear raytrace_3d

        for ii = 1:1:num_elevs
            in_range = IONS.chk_dist(range, ray_N(ii), wgs84);
            
            if debug
                plot3(ray_N(ii).lat,  ray_N(ii).lon, ray_N(ii).height, plotcol, 'HandleVisibility', 'off')
                
                this_nhops = length(ray_data_N(ii).lat);
                for hop = 1:this_nhops
                    if ray_data_N(ii).ray_label(hop) == 1
                        lat = ray_data_N(ii).lat(hop);
                        lon = ray_data_N(ii).lon(hop);
                        plot3(lat, lon, 0, 'x', 'color', plotcol, 'markersize', 20, 'LineWidth',4, 'HandleVisibility', 'off')
                    else
%                         fprintf(ray_data_N(ii).ray_label(hop))
                    end
                end
            else
                if in_range == 1   
                    plot3(ray_N(ii).lat,  ray_N(ii).lon, ray_N(ii).height, plotcol, 'HandleVisibility', 'off')

                    for hop = 1:nhops
                        lat = ray_data_N(ii).lat(hop);
                        lon = ray_data_N(ii).lon(hop);
                        plot3(lat, lon, 0, 'x', 'color', plotcol, 'markersize', 20, 'LineWidth',4, 'HandleVisibility', 'off')
                    end
                elseif in_range == 0
                    fprintf("Elev " + elevs(ii) + " not in range. \n")
                elseif in_range == -1
                    fprintf("Elev " + elevs(ii) + " too high. \n")
                end
            end
        end
    end
    
    view(rot,45)

    % From Python:
    TX = [40.68, -105.04, 0];
    BP = [41.75, -89.62, 200];
    RX = [40.74, -74.18, 0];

    points = [TX; RX];
    fs = 30;
    
    for pt = 1:1:2
        x = points(pt,1); y = points(pt,2); z = points(pt,3);
        txt = "[" + x + ", " + y + ", " + z + "]";
        plot3(x, y, z, 'o', 'color', 'k', 'markersize', 20, 'LineWidth',4, 'HandleVisibility', 'off')
        text(x + 0.25, y, z + 10, txt, 'FontSize', fs-10, 'FontWeight', 'bold') 
    end

    x = BP(1); y = BP(2); z = BP(3);
    txt = "Great Circle Path Midpoint";
    plot3(x, y, z, 'o', 'color', 'k', 'markersize', 20, 'LineWidth',4, 'HandleVisibility', 'off')
    text(x + 0.25, y, z + 10, txt, 'FontSize', fs-10, 'FontWeight', 'bold') 
    
%   Bounds for ionosphere_500km_series_plus
    % grid_corner = [39, -106];
    % plot3(grid_corner(1), grid_corner(2), 0, 'x', 'color', 'k', 'markersize', 20, 'LineWidth',4, 'HandleVisibility', 'off')
    % plot3(grid_corner(1)+4, grid_corner(2), 0, 'x', 'color', 'k', 'markersize', 20, 'LineWidth',4, 'HandleVisibility', 'off')
    % plot3(grid_corner(1), grid_corner(2)+36, 0, 'x', 'color', 'k', 'markersize', 20, 'LineWidth',4, 'HandleVisibility', 'off')
    % plot3(grid_corner(1)+4, grid_corner(2)+36, 0, 'x', 'color', 'k', 'markersize', 20, 'LineWidth',4, 'HandleVisibility', 'off')

    set(gca, 'XDir','reverse')
    hold off
    grid on
    ti = "3D Raytracing Rays from WWV to K2MFF";
    su = UT(i,1) + "-" + UT(i,2) + "-" + UT(i,3) + " | " + ...
         "" + UT(i,4) + ":00Z | " + ...
         "TX: [" + TX(1) + ", " + TX(2) + ...
         "] --> RX: [" + RX(1) + ", " + RX(2) + "] | " + ...
         "Midpoint : [" + BP(1) + ", " + BP(2) + "]";
    
    title([ti,su], 'FontSize', fs)
    xlabel('Latitude (°)', 'FontSize', fs)
    ylabel('Longitude (°)', 'FontSize', fs)
    zlabel('Height (km)', 'FontSize', fs)
    ax = gca;
    ax.FontSize = fs/1.5; 
    
    [i, rot] = update_fig(i, rot);
    
end

function [i, rot] = update_fig(i, rot)

    fig = gcf;
    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'downarrow')
        if i > 1
            i = i - 1;
        end
        clf
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'uparrow')
        i = i + 1;
        clf
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
        rot = rot + 10;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
        rot = rot - 10;
    else
        
    end
    fprintf("\n")
    
    set(gcf, 'Position',  [100, 100, 900, 400])
    
end