clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")

clf

mode_keys = ["O", "X"];
R12_sel = [-1, 25, 50, 57, 100, 200];

count = 1;
for mode_key_i = 1:1:1
    for r12_i = 1:1:1
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
        
        %.ground_range          - geodetic (WGS84) ground range (Km) 
        %.geometric_path_length - Geometrical distance travelled by ray (km)
        
        props = [["ground_range", "ray_data"];
                 ["geometric_path_length", "ray_data"]];
        rps = obj.ray_props(props);

        ground_range = rps.(props(1));
        geometric_path_length = rps.(props(2));

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Plotting
        fprintf("\nPlotting Figure " + count + "\n")
        fig = figure(count);
        tcl = tiledlayout(fig, 2,2);
        h = gobjects(4,1);

        hr_range = 0:1:23;
        elevs_range = flip(elevs);

        nhop_max = 4;
        for nhops = 1:1:nhop_max
            %subplot(2,2,nhops)
            ax = nexttile(tcl);

            hop_field = "hop_" + nhops;
            fprintf("Plotting " + nhops + "-hop rays \n")

            ratio = ground_range.(hop_field) ./ geometric_path_length.(hop_field);
            ratio(isnan(ratio)) = 0;
            ratio = ratio.';
            ratio = flip(ratio);

            %h(nhops) = heatmap(hr_range, elevs_range, ratio, 'ColorLimits', [0.7 1.0], ...
            %            'Colormap', jet);
            
            h(nhops) = heatmap(hr_range, elevs_range, ratio, 'ColorbarVisible', 'off');

            %warning('off', 'MATLAB:structOnObject')
            %hs = struct(h);
            %ylabel(hs.Colorbar, "Ground Range (km) / Geometric Path Length (km)");

            %h(nhops).Title = "Number of Hops:" + nhops;
            %h(nhops).XLabel = 'Time (UT)';
            %h(nhops).YLabel = 'Elevation (°)';
            %h(nhops).GridVisible = 'off';

            % Convert each number in the array into a string
            %CustomYLabels = string(elevs_range);
            % Replace all but the 10th elements by spaces
            %CustomYLabels(mod(elevs_range,10) ~= 0) = " ";
            % Set the 'XDisplayLabels' property of the heatmap 
            % object 'h' to the custom x-axis tick labels
            %h.YDisplayLabels = CustomYLabels;

        end
        
        % Equate color limits in all heatmaps
        globalColorLim = [0, 1];
        set(h, 'ColorLimits', globalColorLim)
        
        % Create global colorbar that uses the global color limits
        ax = axes(tcl,'visible','off','Colormap',h(1).Colormap,'CLim',globalColorLim);
        cb = colorbar(ax);
        %cb.Layout.Tile = 'East';
        cb.Location = 'East';

        ti = "Ground Range / Geometric Path Length Per Elevation Per Hour";
        sgtitle(ti+elevs_string+r12_string+mode_string)
        %set(gcf, 'Position', get(0, 'Screensize') / 1.1);
        
        dirname = "path_length_heatmaps_end/";
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



