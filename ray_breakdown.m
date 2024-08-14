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
            el_inc = 0.2;
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
        
        brk = true;
        obj = IONS(date, elevs, freq, R12, mode, gen, brk);

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Raytracing

        breakdown = obj.get_ray_breakdown();
        per_tot = breakdown.per_tot;
        per_hop = breakdown.per_hop;

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Plotting

        %clf
        hr_range = 0:1:23;
        tmp = zeros(1,5);
        bars = repmat(tmp,24,1);
        nhops = obj.nhops_max;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Aggregate Data
        for hour = 1:1:24
            per_hr = tmp;
            per_hr(1) = per_tot(hour);

            for hop = 1:1:nhops
                hop_field = "hop_" + hop;
                per_hr(hop+1) = per_hop.(hop_field)(hour);
            end

            bars(hour, :) = per_hr;
        end
        
        figure(count)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot Bar Graph
        subplot(1,2,1);

        b = bar(hr_range, bars);
        c = ['k', 'r', 'g', 'b', 'm'];
        for color = 1:1:5
            b(color).FaceColor = c(color);
            b(color).EdgeColor = 'none';
        end

        set(b, {'DisplayName'}, {'Total','1-hop','2-hop','3-hop','4-hop'}');
        legend();
        xlabel('Time (UT)');
        ylabel('Percent of Rays Sent (%)');
        xticks(hr_range);
        ylim([0,0.25])
        grid on;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot Line Graph
        subplot(1,2,2);

        hold on;
        yrange = bars(:, 1).';
        plot(hr_range, yrange, '-k', "LineWidth", 4);
        for hop = 2:1:5
            yrange = bars(:, hop).';
            plot(hr_range, yrange, c(hop), "LineWidth", 2);
        end
        hold off;

        legend({'Total','1-hop','2-hop','3-hop','4-hop'});
        xlabel('Time (UT)');
        ylabel('Percent of Rays Sent (%)');
        xticks(hr_range);
        ylim([0,0.25])
        grid on;

        ti = "Percentage of Recieved Rays by Number of Hops";
        sgtitle(ti+elevs_string+r12_string+mode_string)
        
        count = count + 1;
    end
end
