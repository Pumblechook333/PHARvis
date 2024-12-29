clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")

clf

mode_keys = ["O", "X"];
R12_sel = [-1, 25, 50, 100, 200];

count = 1;
for mode_key_i = 1:1:2
    for r12_i = 1:1:5
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % GET necessary vars

        date = [2021 7 1 0 0];

        el_start = 0; 
        
        hi_res = 0;
        if hi_res
            el_inc = 0.1;
            el_stop = 50;
        else
            el_inc = 5;
            el_stop = 90;
        end
        
        elevs = el_start : el_inc : el_stop;

        freq = 10;
        R12 = R12_sel(r12_i);
        mode_key = mode_keys(mode_key_i);
        mode_map = struct('O', 1, 'No', 0, 'X', -1);
        mode = mode_map.(mode_key);
        gen = 0; % 0 = no gen, 1 = gen

        elevs_string = " || Initial Elevations: " ...
                       + el_start + ":" + el_inc + ":" + el_stop;
        r12_string = " || R12: " + R12;
        mode_string = " || " + mode_key + "-mode";
        
        obj = IONS(date, elevs, freq, R12, mode, gen);

        iono = obj.get_iono_parms();

        g = obj.get_gen_params();
        elevs = g.elevs;

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Raytracing

        %.phase_path              - phase path (km) \
        %                         ** NUMBER OF PHASES THE SIGNAL GONE THROUGH
        %                           (1 cycle, 2 cycles, 2.75 cycles)
        %.geometric_distance      - geometrical distance travelled by ray (km)
        
        props = [["phase_path", "ray"];
                 ["geometric_distance", "ray"]];
        rps = obj.ray_props(props);

        phase_path = rps.(props(1));
        geometric_distance = rps.(props(2));

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Plotting
        fprintf("\nPlotting Figure: " + count + "\n")
        
        figure(count)

        hr_range = 0:1:23;
        elevs_range = flip(elevs);

        nhop_max = 4;
        for nhops = 1:1:nhop_max
            subplot(2,2,nhops)

            hop_field = "hop_" + nhops;
            fprintf("Plotting " + nhops + "-hop rays \n")

            ratio = phase_path.(hop_field) ./ geometric_distance.(hop_field);
            ratio(isnan(ratio)) = 0;
            ratio = ratio.';
            ratio = flip(ratio);

            h = heatmap(hr_range, elevs_range, ratio, 'ColorLimits', [0.7 1.0], ...
                        'Colormap', jet);

            warning('off', 'MATLAB:structOnObject')
            hs = struct(h);
            ylabel(hs.Colorbar, "Phase Path (m) / Geometric Distance (m)");

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
        sgtitle(ti+elevs_string+r12_string+mode_string)
        
        set(gcf, 'Position', get(0, 'Screensize') / 1.1);
        
        dirname = "phase_path_heatmaps/";
        if not(isfolder(dirname))
                    mkdir(dirname)
        end
        
        if hi_res
            set(gcf,'visible','off')
            figname = "figure_" + R12 + "_" + mode_key + ".jpg";
            sppi = get(groot,"ScreenPixelsPerInch");
            exportgraphics(gcf, dirname+figname, 'Resolution', sppi)
        end

        count = count +1;
    end
end
