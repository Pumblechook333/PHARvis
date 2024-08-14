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

hi_res = 1;
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
brk = true;

%% Loop R12 densities
count = 1;
r12_max = 5;
for r12_i = 1:1:r12_max
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET necessary vars
    R12 = R12_sel(r12_i);

    elevs_string = " || Initial Elevations: " ...
                   + el_start + ":" + el_inc + ":" + el_stop;
    r12_string = " || R12: " + R12;

    obj_O = IONS(date, elevs, freq, R12, 1, gen, brk);
    obj_X = IONS(date, elevs, freq, R12, -1, gen, brk);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Raytracing

    breakdown_O = obj_O.get_ray_breakdown();
    per_tot_O = breakdown_O.per_tot;
    per_hop_O = breakdown_O.per_hop;
    
    breakdown_X = obj_X.get_ray_breakdown();
    per_tot_X = breakdown_X.per_tot;
    per_hop_X = breakdown_X.per_hop;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting

    %clf
    hr_range = 0:1:23;
    tmp = zeros(1,5);
    bars_O = repmat(tmp,24,1);
    bars_X = repmat(tmp,24,1);
    nhops = obj_O.nhops_max;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aggregate Data
    for hour = 1:1:24
        per_hr_O = tmp;
        per_hr_X = tmp;
        
        per_hr_O(1) = per_tot_O(hour);
        per_hr_X(1) = per_tot_X(hour);

        for hop = 1:1:nhops
            hop_field = "hop_" + hop;
            per_hr_O(hop+1) = per_hop_O.(hop_field)(hour);
            per_hr_X(hop+1) = per_hop_X.(hop_field)(hour);
        end

        bars_O(hour, :) = per_hr_O;
        bars_X(hour, :) = per_hr_X;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot Line Graph
    
    figure(count)
    hold on;
    yrange_O = bars_O(:, 1).';
    yrange_X = bars_X(:, 1).';
    
    c = ["k", "r", "g", "b", "m"];
    
    plot(hr_range, yrange_O, '-k', "LineWidth", 4);
    plot(hr_range, yrange_X, '--k', "LineWidth", 4);
    for hop = 2:1:5
        style = "-" + c(hop);
        yrange_O = bars_O(:, hop).';
        plot(hr_range, yrange_O, style, "LineWidth", 2);
        
        style = "--" + c(hop);
        yrange_X = bars_X(:, hop).';
        plot(hr_range, yrange_X, style, "LineWidth", 2);
    end
    hold off;
    
    legend_cells = {'Total O', 'Total X', '1-hop O', '1-hop X',...
                    '2-hop O', '2-hop X','3-hop O', '3-hop X',...
                    '4-hop O', '4-hop X'};
    legend(legend_cells, 'Location', 'eastoutside');
    xlabel('Time (UT)');
    ylabel('Percent of Rays Sent (%)');
    xticks(hr_range);
    ylim([0,0.25])
    grid on;
    
    set(gca,"FontSize",20)

    ti = "Percentage of Recieved Rays by Number of Hops";
    title(ti+elevs_string+r12_string)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE PLOTS if high res
    set(gcf, 'Position', get(0, 'Screensize') / 1.1);
    
    dirname = "breakdown_plots_OX/";
    if not(isfolder(dirname))
                mkdir(dirname)
    end
    
    if hi_res
        set(gcf,'visible','off')
        figname = "figure_" + R12 + "_OX" + ".jpg";
        sppi = get(groot,"ScreenPixelsPerInch");
        exportgraphics(gcf, dirname+figname, 'Resolution', sppi)
    end

    count = count + 1;
end
