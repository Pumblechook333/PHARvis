clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n")

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% GET necessary vars

date = [2021 7 1 0 0];

el_start = 0;
el_inc = 0.2;
el_stop = 50;
% el_inc = 5;
% el_stop = 90;
elevs = el_start : el_inc : el_stop;

freq = 10;
R12_sel = [-1, 25, 50, 100, 200];
R12 = R12_sel(4);
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

props = [["ground_range", "ray_data"];
         ["geometric_path_length", "ray_data"]];
rps = obj.ray_props(props);

ground_range = rps.(props(1));
geometric_path_length = rps.(props(2));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plotting

clf

hr_range = 0:1:23;
elevs_range = flip(elevs);

nhop_max = 4;
for nhops = 1:1:nhop_max
    subplot(2,2,nhops)

    hop_field = "hop_" + nhops;
    fprintf("Plotting " + nhops + "-hop rays \n")
    
    ratio = ground_range.(hop_field) ./ geometric_path_length.(hop_field);
    ratio(isnan(ratio)) = 0;
    ratio = ratio.';
    ratio = flip(ratio);
    
    h = heatmap(hr_range, elevs_range, ratio, 'ColorLimits', [0.7 1.0], ...
                'Colormap', jet);
            
    warning('off', 'MATLAB:structOnObject')
    hs = struct(h);
    ylabel(hs.Colorbar, "Ground Range (km) / Geometric Path Length (km)");
    
    h.Title = "Number of Hops:" + nhops;
    h.XLabel = 'Time (UT)';
    h.YLabel = 'Elevation (°)';
    h.GridVisible = 'off';
    
    % Convert each number in the array into a string
    CustomYLabels = string(elevs_range);
    % Replace all but the 10th elements by spaces
    CustomYLabels(mod(elevs_range,10) ~= 0) = " ";
    % Set the 'XDisplayLabels' property of the heatmap 
    % object 'h' to the custom x-axis tick labels
    h.YDisplayLabels = CustomYLabels;

end

ti = "Ground Range / Geometric Path Length Per Elevation Per Hour";
sgtitle(ti+elevs_string+r12_string)



