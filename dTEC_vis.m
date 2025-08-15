% clc;
% clear;

datfile = "E:\Sabastian\Perry_Lab\dTEC_data\tid_182_2021_w121_n01_e15.h5";
tecs_size = [10, 161926918];
if ~exist('filtered') filtered = struct(); end

%%

h5disp(datfile, '/tecs')

%%

% y elements per row, y rows
x = 10;
y = 100;
st = [1,1];
en = [x,y];

if ~exist('ak')
    ak = h5read(datfile,'/tecs', st, en); 
    ak = ak(:,2:end)';
end

%%

doyut=ak(:,1); 
ut=(doyut-fix(doyut)).*24;
doy = fix(doyut) + 1;           % day number of the year: note doyut starts 0.0 as day 1 at 0 UT
la=ak(:,2);                     % pierce point latitude
lo=ak(:,3);                     % pierce point longitude
tid=ak(:,4);                    % dtec
tec=ak(:,5);                    % tec  
egf=ak(:,6);                    % edge flag
prn=ak(:,7);                    % PRN
sit=ak(:,8);                    % SITE
ele=ak(:,9);                    % elevation
sat=ak(:,10);                   % GPS=0;  GLONASS=1

%%

chunksize = 1000000;      % Select how much of column to store

ncycles = round(tecs_size(2)/chunksize) - 1;

% Check 2, retrieve   1,  3,  and 4
% Check lat, retrieve ut, lon and dTEC

for i = 1:1:ncycles

    ch_name = "ch_" + num2str(i);
    fnames = fieldnames(filtered);
    
    if ~ismember(ch_name,fnames,'rows')
        fprintf("Chunk " + i + " / " + ncycles + "\n")
        chunk = struct();
    
        lats = sel_hdf5_chunk(datfile, i, chunksize, 2);
        
        uts = sel_hdf5_chunk(datfile, i, chunksize, 1);
        lons = sel_hdf5_chunk(datfile, i, chunksize, 3);
        dtecs = sel_hdf5_chunk(datfile, i, chunksize, 4);
    
        if ~(i == 1) len = chunksize; else len = (chunksize - 1); end
        
        filt_lats =     [];
        filt_lons =     [];
        filt_uts =      [];
        filt_dtecs =    [];
        for j = 1:1:len
            if (lats(j) > 40) && (lats(j) < 42)
                filt_lats =     [filt_lats lats(j)];
                filt_lons =     [filt_lons lons(j)];
                filt_uts =      [filt_uts  uts(j)];
                filt_dtecs =    [filt_dtecs dtecs(j)];
            end
        end
        
        chunk.lats = filt_lats';
        chunk.lons = filt_lons';
        chunk.uts = filt_uts';
        chunk.dtecs = filt_dtecs';
        
        filtered.(ch_name) = chunk;
    end

end


%%

filtered_arr = [];

comb_lats = [];
comb_lons = [];
comb_uts = [];
comb_dtecs = [];

for i = 1:1:ncycles
    ch_name = "ch_" + num2str(i);
    fprintf("Chunk " + i + " / " + ncycles + "\n")
    
    comb_lats = [comb_lats; filtered.(ch_name).lats];
    comb_lons = [comb_lons; filtered.(ch_name).lons];
    comb_uts = [comb_uts; filtered.(ch_name).uts];
    comb_dtecs = [comb_dtecs; filtered.(ch_name).dtecs];

end

filtered_arr = [comb_uts, comb_lats, comb_lons, comb_dtecs];

writematrix(filtered_arr, 'filtered_dtec.csv');

%%


function selection = sel_hdf5_chunk(filename, iter, chunksize, col)

    st = [col,      1 + chunksize*(iter-1)];
    en = [1, chunksize];
    
    selection = h5read(filename,'/tecs', st, en); 
    if iter == 1
        selection = selection(:,2:end)';
    else
        selection = selection(:,1:end)';
    end

end


