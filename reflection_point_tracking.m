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
        
        hi_res = 1;
        if hi_res
            el_inc = 0.1;
            el_stop = 50;
        else
            el_inc = 5;
            el_stop = 90;
        end
        
        elevs = el_start : el_inc : el_stop;
        num_elevs = length(elevs);

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
        iono_height = iono.iono_height;

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Raytracing

        props = [["height", "ray_max"]];
        rps = obj.ray_props(props);
        max_heights = rps.(props(1));

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Plotting
        fprintf("\nPlotting Figure: " + count + "\n")
        
        clf
        figure(count)
        set(gcf,'visible','off')
        
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

        for nhops = 1:1:obj.nhops_max
            hop_field = "hop_" + nhops;
            fprintf("Plotting " + nhops + "-hop rays \n")

            hr_range = 0:1:24;

            for k = 1:1:24
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
        title((ti+elevs_string+r12_string+mode_string), 'FontSize', fs)
        xlabel('Time (UT)', 'FontSize', fs)
        ylabel('Height (km)', 'FontSize', fs)

        xticks(0:1:24)

        ax = gca;
        ax.FontSize = fs/1.5;
        
        set(gcf, 'Position', get(0, 'Screensize') / 1.1);
                
        dirname = "reflection_point_plots/";
        if not(isfolder(dirname))
                    mkdir(dirname)
        end
        
        if hi_res
            figname = "figure_" + R12 + "_" + mode_key + ".jpg";
            sppi = get(groot,"ScreenPixelsPerInch");
            exportgraphics(gca, dirname+figname, 'Resolution', sppi)
        end
        count = count + 1;
    end
end
