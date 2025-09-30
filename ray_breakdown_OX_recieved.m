%% Setup
clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")

clf
%% Constants / Settings
mode_keys = ["O", "X"];
% r12_sel = [-1, 25, 50, 100, 200];
% r12_sel = [132, 145, 137, 113, 111, 104, 97, 74, 66, 57, 67]; % 2023 10 10-20
% r12_sel = [78, 80, 84, 77, 55, 55, 90, 86, 120, 148, 171]; % 2024 04 05-15
% r12_sel = [78, 80];
r12_sel = 57;

r12_sz = size(r12_sel);
r12_sz = r12_sz(2);

year = 2021;
month = 7;
days = [1];

% year = 2023;
% month = 10;
% days = 10:1:20;

% year = 2024;
% month = 04;
% days = 5:1:15;
% days = [5, 6];

el_start = 0;

hi_res = 0;
if hi_res
    el_inc = 0.2;
    el_stop = 50;
else
    el_inc = 2;
    el_stop = 50;
end

elevs = el_start : el_inc : el_stop;

freq = 10;

gen = 0; % 0 = no gen, 1 = gen
brk = true;

%% Loop R12 densities
count = 1;
r12_max = r12_sz;
for r12_i = 1:1:r12_max
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET necessary vars
    R12 = r12_sel(r12_i);
    date = [year month days(r12_i) 0 0];
    
    date_string = " || " + year + "-" + month + "-" + days(r12_i);
    r12_string = " || R12: " + R12;
    elevs_string = " || Initial Elevations: " ...
                   + el_start + ":" + el_inc + ":" + el_stop;
    
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
    hr_range = 0:1:24;
    tmp = zeros(1,4);           % One slot for each hop
    bars_O = repmat(tmp,25,1);
    bars_X = repmat(tmp,25,1);
    nhops = obj_O.nhops_max;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aggregate Data
    for hour = 1:1:25
        fprintf("Aggregating Data for hour: " + hour + "\n")
        per_hr_O = tmp;
        per_hr_X = tmp;

        for hop = 1:1:nhops
            hop_field = "hop_" + hop;
            
            if per_tot_O(hour) == 0
                per_hr_O(hop) = 0;
            else
                per_hr_O(hop) = per_hop_O.(hop_field)(hour) / per_tot_O(hour);
            end
            
            if per_tot_X(hour) == 0
                per_hr_X(hop) = 0;
            else
                per_hr_X(hop) = per_hop_X.(hop_field)(hour) / per_tot_X(hour);
            end
            
        end

        bars_O(hour, :) = per_hr_O;
        bars_X(hour, :) = per_hr_X;
    end
    
    if hi_res
        fprintf("Exporting data for R12: " + R12 + "\n\n")
        writematrix(bars_O, "EXPORT/export_data/" + R12 + "_O_mode_percentages.csv")
        writematrix(bars_X, "EXPORT/export_data/" + R12 + "_X_mode_percentages.csv")
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot Line Graph
    
    figure(count)
    hold on;
    
    c = ["r", "g", "b", "m"];
    
    for hop = 1:1:4
        style = "-" + c(hop);
        yrange_O = bars_O(:, hop).';
        plot(hr_range, yrange_O, style, "LineWidth", 2);
        
        style = "--" + c(hop);
        yrange_X = bars_X(:, hop).';
        plot(hr_range, yrange_X, style, "LineWidth", 2);
    end
    hold off;
    
    legend_cells = {'1-hop O', '1-hop X','2-hop O', '2-hop X',...
                    '3-hop O', '3-hop X','4-hop O', '4-hop X'};
    legend(legend_cells, 'Location', 'eastoutside');
    xlabel('Time (UT)');
    ylabel('Percent of Rays Recieved (%)');
    xticks(hr_range);
    xlim([0,24]);
    ylim([0,1])
    grid on;
    
    set(gca,"FontSize",20)

    ti = "Ray Hop Breakdown";
    title(ti+date_string+elevs_string+r12_string)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE PLOTS if high res
    set(gcf, 'Position', get(0, 'Screensize') / 1.1);
    
    dirname = "PLOTS/breakdown_plots_OX_recieved/";
    if not(isfolder(dirname))
        mkdir(dirname)
    end
    
    if hi_res
        set(gcf,'visible','off')
        sppi = get(groot,"ScreenPixelsPerInch");
        figname = "figure_" + R12 + "_OX_recieved" + ".jpg";
        exportgraphics(gcf, dirname+figname, 'Resolution', sppi)
    end

    count = count + 1;
end
