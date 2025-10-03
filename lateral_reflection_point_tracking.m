% Exports .csv files containing the following information for for 
% successfully recieved rays:
% - maximum height of simulated ray
% - latitude and longitude coordinates for where maximum height is acheived

%%
% Clear console output and memory

clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")
clf

%%
% Performs simulation for O and X modes at provided R12 value(s)

mode_keys = ["O", "X"];
r12_sel = [57];
r12_sz = size(r12_sel);
r12_sz = r12_sz(2);
count = 1;

r12_max = r12_sz;
for mode_key_i = 1:1:2
    for r12_i = 1:1:r12_max
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
            el_stop = 50;
        end
        
        elevs = el_start : el_inc : el_stop;
        num_elevs = length(elevs);

        freq = 10;
        R12 = r12_sel(r12_i);
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
        
        props = [["height", "ray_max"]; ["lat", "ray_data_all"]; 
                 ["lon", "ray_data_all"]];
        rps = obj.ray_props(props);
        max_heights = rps.(props(1));
        lats = rps.(props(2));
        lons = rps.(props(3));
        
        for nhops = 1:1:obj.nhops_max       % Per Hop
            hop_field = "hop_" + nhops;
            fprintf("Plotting " + nhops + "-hop rays \n")
            
            mh = max_heights.(hop_field);
            la = lats.(hop_field);
            lo = lons.(hop_field);
            
            mh_str = cellArrayToStringMatrix(mh);
            la_str = cellArrayToStringMatrix(la);
            lo_str = cellArrayToStringMatrix(lo);
            
            if hi_res
                fprintf("Exporting data for R12: " + R12 + " / " + mode_key + "-mode / " + nhops + "-hop \n\n")
                writecell(mh_str, "EXPORT/max_heights_" + R12 + "_" + mode_key + "mode_" + nhops + "nhops.csv")
                writecell(la_str, "EXPORT/lats_" + R12 + "_" + mode_key + "mode_" + nhops + "nhops.csv")
                writecell(lo_str, "EXPORT/lons_" + R12 + "_" + mode_key + "mode_" + nhops + "nhops.csv")
            end
            
        end
                
        count = count + 1;
    end
end


%% 
% Function to convert cell array output from IONS into format exportable 
% to .csv format 

function stringMatrix = cellArrayToStringMatrix(cellArray)

if ~iscell(cellArray)
    warning('Input is not a cell array. Returning empty string.');
    stringMatrix = char.empty(); % Return empty char array
    return;
end

if isempty(cellArray)
    stringMatrix = char.empty();
    return;
end

[rows, cols] = size(cellArray);  % Get dimensions of the cell array
stringMatrix = cell(rows, cols); % Preallocate the cell array of strings

for i = 1:rows
    for j = 1:cols
        currentValue = cellArray{i, j};

        if ischar(currentValue)
            stringMatrix{i, j} = currentValue;
        elseif isnumeric(currentValue)
            stringMatrix{i, j} = num2str(currentValue);
        elseif islogical(currentValue)
            stringMatrix{i, j} = num2str(currentValue);
        else
            warning(['Unsupported data type at (' num2str(i) ', ' num2str(j) '). Converting to string ''Unsupported''.']);
            stringMatrix{i, j} = 'Unsupported';  % Or handle differently
        end
    end
end

end
