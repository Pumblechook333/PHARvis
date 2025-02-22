clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")

clf

mode_keys = ["O", "X"];
% R12_sel = [-1, 25, 50, 57, 100, 200];
R12_sel = [57];

count = 1;
for mode_key_i = 1:1:2
    for r12_i = 1:1:1
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % GET necessary vars

        date = [2021 7 1 0 0];

        el_start = 0;
        
        hi_res = 1;
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
        tcl.Padding = 'compact';
        tcl.TileSpacing = 'compact';
        ax = gobjects(4,1); % Store axes handles
        h = gobjects(4,1);

        hr_range = 0:1:24;
        elevs_range = flip(elevs);

        nhop_max = 4;
        for nhops = 1:1:nhop_max
            ax(nhops) = nexttile(tcl); % Get axes handle and store it

            hop_field = "hop_" + nhops;
            fprintf("Plotting " + nhops + "-hop rays \n")

            gr = cell2mat(ground_range.(hop_field));
            pl = cell2mat(geometric_path_length.(hop_field));

            ratio = gr ./ pl;
            ratio(isnan(ratio)) = 0;
            ratio = ratio.';
            ratio = flip(ratio);

            h(nhops) = imagesc(hr_range, elevs_range, ratio);
            axis xy;
            colormap(ax(nhops), jet); % Set colormap for each axes

            title(ax(nhops), "Number of Hops:" + nhops); % Set title for the axes
            xlabel(ax(nhops), 'Time (UT)');
            ylabel(ax(nhops), 'Elevation (°)');

            % Custom Y-axis labels (improved)
            CustomYLabels = string(elevs_range);
            CustomYLabels(mod(elevs_range,10) ~= 0) = " ";
            ax(nhops).YTick = flip(elevs_range(mod(elevs_range,10) == 0));
            ax(nhops).YTickLabel = flip(CustomYLabels(mod(elevs_range,10) == 0));
            
             % Custom X-axis labels (improved)
            CustomXLabels = string(hr_range);
            CustomXLabels(mod(hr_range,2) ~= 0) = " ";
            ax(nhops).XTick = hr_range(mod(hr_range,2) == 0);
            ax(nhops).XTickLabel = CustomXLabels(mod(hr_range,2) == 0);
            
            ax(nhops).FontSize = 14;

        end
        
        % Equate color limits in all heatmaps (Corrected)
%         data = cell(1, nhop_max);
%         for nhops = 1:nhop_max
%             data{nhops} = get(h(nhops), 'CData');
%         end
% 
%         globalColorLim = [min(cellfun(@(x) min(x(:)), data)), max(cellfun(@(x) max(x(:)), data))];
        globalColorLim = [0.7, 1.0];
        
        for nhops = 1:nhop_max  % Loop through axes handles
            set(ax(nhops), 'CLim', globalColorLim); % Set CLim for each axes
        end

        % Create global colorbar to the right of the tiled layout
        cb_ax = axes(fig, 'Position', [0.96 0.1 0.03 0.8]);
        set(cb_ax, 'XTickLabel', {}, 'YTickLabel', {})
        set(cb_ax,'color','none')
        cb_ax.XAxis.Visible = 'off';
        cb_ax.YAxis.Visible = 'off';
        
        colormap(cb_ax, jet);
        caxis(globalColorLim);
        cb = colorbar(cb_ax);
        set(cb, 'Position', [0.96 0.1 0.02 0.8])
        set(cb, 'YAxisLocation','right', 'FontSize', 10)

        % Improve Title
        ti = "Ground Range / Geometric Path Length Per Elevation Per Hour";
        sgtitle(tcl, ti+elevs_string+r12_string+mode_string, 'Fontsize', 18)

        dirname = "path_length_heatmaps_end/";
        if not(isfolder(dirname))
                    mkdir(dirname)
        end
        
        set(gcf, 'Position',  [100, 100, 1400, 700])
        
        if hi_res
            set(gcf,'visible','off')
            figname = "figure_" + R12 + "_" + mode_key + ".jpg";
            sppi = get(groot,"ScreenPixelsPerInch");
            exportgraphics(gcf, dirname+figname, 'Resolution', sppi)
        end

        count = count +1;
    end
end



