clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")

clf

mode_keys = ["O", "X"];
% r12_sel = [57];
% r12_sz = size(r12_sel);
% r12_sz = r12_sz(2);
count = 1;

% 2023 10 10-20
% year = 2023;
% month = 10;
% days = 10:1:20;
% R12s = [132, 145, 137, 113, 111, 104, 97, 74, 66, 57, 67]; % 2023 10 10-20

% 2024 04 05-15
year = 2024;
month = 04;
days = 5:1:15;
R12s = [78, 80, 84, 77, 55, 55, 90, 86, 120, 148, 171];

% r12_max = r12_sz;
for mode_key_i = 1:1:2
    % for day = 1:1:10
    for day = 5:1:5
        % for r12_i = 1:1:r12_max
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % GET necessary vars
    
        date = [year month days(day) 0 0];
    
        el_start = 0;
        
        hi_res = 0;
        if hi_res
            el_inc = 0.1;
            el_stop = 50;
        else
            el_inc = 2;
            el_stop = 90;
        end
        
        elevs = el_start : el_inc : el_stop;
        num_elevs = length(elevs);
    
        freq = 10;
        % R12 = r12_sel(r12_i);
        R12 = R12s(day);
        % R12 = 57;
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
    
        props = [["height", "ray_max"]; ["initial_elev", "ray"]];
        rps = obj.ray_props(props);
        max_heights = rps.(props(1));
        initial_elevs = rps.(props(2));
    
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Plotting
        fprintf("\nPlotting Figure: " + count + "\n")
        
    %         clf
        figure(count)
    %         set(gcf,'visible','off')
        
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
    
        for nhops = 1:1:obj.nhops_max       % Per Hop
            hop_field = "hop_" + nhops;
            fprintf("Plotting " + nhops + "-hop rays \n\n")
            
            for k = 1:1:25                  % Per Hour
                hour_field = 'i' + string(k);
    
                for i = 1:1:num_elevs       % Per Elevation
                    y = max_heights.(hop_field){k,i};
                    % in_range = (0 < y) & (y < iono_height)
                    initial_elev = initial_elevs.(hop_field){k,i};

                    in_range = ~(y == 0);
                    
                    if in_range
                        sz = initial_elev / el_stop;    % Size the markers with respect to initial elevation
                        if sz == 0, sz = 1; end

                        scaled_size = 100*sz;
                        plot(k-1, y, 'color', colors(nhops), 'marker', '.', ...
                            'markersize', scaled_size, 'LineWidth', 4)
                    end
                end
    
            end
            
            if hi_res
                fprintf("Exporting data for R12: " + R12 + " / " + mode_key + "-mode / " + nhops + "-hop \n\n")
                writematrix(max_heights.(hop_field), "export_data/max_heights_" + R12 + "_" + mode_key + "mode_" + nhops + "nhops.csv")
                writematrix(initial_elevs.(hop_field), "export_data/initial_elevs_" + R12 + "_" + mode_key + "mode_" + nhops + "nhops.csv")
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
                
        dirname = "PLOTS/reflection_point_plots/";
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
