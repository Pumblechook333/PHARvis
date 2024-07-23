clear
clc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Ionospheric Grid Parameters

% UT = [2021 7 1 10 32];          % Sunrise
UT = [2021 7 1 1 34];           % Sunset
R12 = -1;
doppler_flag = 1;               % interested in Doppler shift

% Bottom left corner of grid (lat, lon) in [deg]
% grid_corner = [40, -(105+(48/60))];
% grid_corner = [38, -106];
grid_corner = [36, -110];

ht_start = 0;          % start height for ionospheric grid (km)
ht_inc = 4;             % height increment (km)
num_ht = 100 +1;     

lat_start = grid_corner(1);
lat_inc = 0.5;
num_lat = 18 +1;

lon_start= grid_corner(2);
lon_inc = 1.0;
num_lon = 100 +1;

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Magnetospheric Grid Parameters

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = ht_inc;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);

B_lat_start = lat_start;
B_lat_inc = lat_inc;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);

B_lon_start = lon_start;
B_lon_inc = lon_inc;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc); 

geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate Ionosphere
fname_append = "_400km.mat";
% fname = 'ionosphere_800km.mat';
fname = 'ionosphere' + fname_append;

gen = 0;
if gen == 1
    fprintf('Generating ionospheric and geomag grids... \n')
    
    tic
    
    [iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, doppler_flag);
    
    save(fname, 'iono_pf_grid', 'iono_pf_grid_5' ,...
        'collision_freq', 'Bx', 'By', 'Bz')
    
    toc
else
    fprintf('Loading ionospheric and geomag grids... \n')
    load = matfile(fname);
    iono_pf_grid = load.iono_pf_grid;
    iono_pf_grid_5 = load.iono_pf_grid_5;
    collision_freq = load.collision_freq;
    Bx = load.Bx;
    By = load.By;
    Bz = load.Bz;
end

fprintf('\n')

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Honing Parameters
tic

TX_coord = [40.67583063, -105.038933178];   % WWV Station
RX_coord = [40.742018, -74.178975];         % K2MFF Station

height = 0; % Assuming from ground
reciever_location = struct('lat', RX_coord(1), 'lon', RX_coord(2), 'alt', height);

% 1e-12 to 1e-2, .001 to 1, 1 to 100
% solver tolerance, minimum step size, maximum step size
% tolerance: distance allowed from step size
% step size: distance between points on ray
% tolerance = [1e-7 0.01 25];
% tolerance = [1e-12 0.005 1];
tolerance = [1e-6 0.5 1]; %%

% Mode to consider:
% 1: 1-hop transmission with reflection at 250 km altitude (F-layer)
% 2: 1-hop transmission with reflection at 120 km altitude (E-layer)
% 3: 2-hop transmission with reflection at 250 km altitude (F-layer)
% 4: 2-hop transmission with reflection at 120 km altitude (E-layer)
% 5: 1-hop transmission with reflection at 80 km altitude  (D-layer)
% 6: 2-hop transmission with reflection at 80 km altitude  (D-layer)
% 7: 4-hop transmission with reflection at 250 km altitude (F-layer)
% 8: 4-hop transmission with reflection at 120 km altitude (E-layer)
% 9: 4-hop transmission with reflection at 80 km altitude  (D-layer)
mode = 9;
% NOTE: Modes 2 and 5 determine negative elevation

% 0-1, 1-100
ide = 0.5; idb = 1; 
% ide = 1; idb = 10; %%

freq = 10;  % Frequency of transmission (10MHz)
OX = 0;     % -1 = Extraordinary, 0 = No field, 1 = Ordinary

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Begin Honing
hone = 1;
if hone==1
    fname_append = fname_append + "_sunset.mat";

    fprintf('Honing rays ...\n\n');

    rece_struct = struct('lat',[],'lon',[],'height',[],'hour',[],'minute',[],...
                         'elevs',[],'bearing',[],'ground_range',[],...
                         'group_range',[],'absorption',[], 'electron_density',[],...
                         'phase_path',[], 'refractive_index',[], 'range',[]);
    data_struct = struct('lat',[],'lon',[],'ground_range',[],...
                         'group_range',[],'phase_path',[],'initial_elev',[],...
                         'final_elev',[],'initial_bearing',[],'final_bearing',[],...
                         'total_absorption',[],'deviative_absorption',[],...
                         'TEC_path',[],'Doppler_shift',[],'apogee',[],...
                         'geometric_path_length',[],'frequency',[],...
                         'nhops_attempted',[], 'ray_label',[],...
                         'NRT_elapsed_time',[]);
    path_struct = struct('initial_elev',[], 'initial_bearing',[],...
                         'frequency',[], 'lat',[], 'lon',[], 'height',[],...
                         'group_range',[], 'phase_path',[],...
                         'refractive_index',[], 'group_refractive_index',[],...
                         'wavenorm_ray_angle',[], 'wavenorm_B_angle',[],...
                         'polariz_mag',[], 'wave_Efield_tilt',[],...
                         'volume_polariz_tilt',[], 'electron_density',[],...
                         'geomag_x',[], 'geomag_y',[], 'geomag_z',[],...
                         'geometric_distance',[], 'collision_frequency',[],...
                         'absorption',[]);
    rays_received = repmat(rece_struct,1,18);
    rays_data = repmat(data_struct,1,18);
    rays_path = repmat(path_struct, 1,18);

    count = 1;
    for rece=1:2
        if rece==1
            rece_dis=1e6;
        end
        if rece==2
            rece_dis=1e5;
        end

        % Rece_dist : maximum distance (m) allowed for landing coordinate from reciever
        % Newton_params = struct('rece_dis', 3e6, 'iteration_num', 10); 
        Newton_params = struct('rece_dis', rece_dis, 'iteration_num', 10); %%

        fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        fprintf("Allowed land distance from receiver: " + rece_dis/1e3 + " km \n\n");
        for mode=1:9
            raytrace_params = struct('origin_lat', TX_coord(1), 'origin_lon', TX_coord(2), 'origin_alt', height,... 
                             'init_delta_elev', ide, 'init_delta_bear', idb, 'fr', freq,...
                             'mode', mode, 'OX_mode', OX, 'tol', tolerance);

            ray_receive = Honing_ground_to_ground(reciever_location, Newton_params,... 
            raytrace_params, iono_grid_parms, geomag_grid_parms, iono_en_grid, ...
            iono_en_grid_5, collision_freq, Bx, By, Bz);

            rays_data(count) = ray_data(1);
            rays_path(count) = ray_path(1);
            rays_received(count) = ray_receive;
            count = count + 1;
        end
    end

    NRT_total_time = toc;
    fprintf('\nTotal mex execution time = %f \n\n', NRT_total_time)

    save("rays_received"+fname_append, "rays_received");
    save("rays_path"+fname_append, "rays_path");
    save("rays_data"+fname_append, "rays_data");
else
    fprintf('Loading rays ...\n\n');
    rays_received = matfile("rays_received"+fname_append);
    rays_received = rays_received.rays_received;
    rays_path = matfile("rays_path"+fname_append);
    rays_path = rays_path.rays_path;
    rays_data = matfile("rays_data"+fname_append);
    rays_data = rays_data.rays_data;
end
    
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plotting

for raynum=2
    fprintf("Ray #" + raynum + ": \n")
    if length(rays_received(raynum).lat) > 1
        fprintf('Plotting received ray... \n\n')

        figure(raynum)
        pos = get(gcf, 'position');
        box_scale = 1.5;
        pos(3) = pos(3)*box_scale;
        pos(4) = pos(4)*box_scale;
        set(gcf, 'position', pos)
        plot3(rays_received(raynum).lat, mod(rays_received(raynum).lon, 360), rays_received(raynum).height, '.b', ...
              'markersize', 5)
        plot3(rays_received(raynum).lat, rays_received(raynum).lon, rays_received(raynum).height, '.b', ...
              'markersize', 5)
        set(gca, 'Zlim', [0 500])
        grid on
        t = sprintf('Honed ray | UTC: %d-%d-%d | From [%f, %f] to [%f, %f]', ...
            UT(1), UT(2), UT(3), TX_coord(1), TX_coord(2), RX_coord(1), RX_coord(2));
        title(t)
        xlabel('latitude (deg)')
        ylabel('longitude (deg)')
        zlabel('Height (km)')
    else
        fprintf('No ray recieved - change initial parameters. \n')

        lbl = rays_data(raynum).ray_label;
        fprintf("Failure: ")
        switch lbl
            case 1
                fprintf("Received")
            case -1
                fprintf("Backscatter")
            case -2
                fprintf("Penetrated ionosphere")
            case -3
                fprintf("Ray exited ionospheric grid")
            case -4
                fprintf("Ray exceeded maximum allowed number of points")
            case -100
                fprintf("CATASTROPHIC ERROR")
        end
        fprintf("\n")

        fprintf('Plotting last traced ray... \n\n')

        figure(raynum)
        pos = get(gcf, 'position');
        box_scale = 1.5;
        pos(3) = pos(3)*box_scale;
        pos(4) = pos(4)*box_scale;
        set(gcf, 'position', pos)
        plot3(rays_path(raynum).lat, mod(rays_path(raynum).lon, 360), rays_path(raynum).height, '.b', ...
              'markersize', 5)
        plot3(rays_path(raynum).lat, rays_path(raynum).lon, rays_path(raynum).height, '.b', ...
              'markersize', 5)
        set(gca, 'Zlim', [0 500])
        grid on
        t = sprintf('Honed ray | UTC: %d-%d-%d | From [%f, %f] to [%f, %f]', ...
            UT(1), UT(2), UT(3), TX_coord(1), TX_coord(2), RX_coord(1), RX_coord(2));
        title(t)
        xlabel('latitude (deg)')
        ylabel('longitude (deg)')
        zlabel('Height (km)')
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function ray_receive=Honing_ground_to_ground(receiver_location,...
    Newton_params,raytrace_params,iono_grid_parms,geomag_grid_parms,Ne,Ne_5,...
    collision_freq,Bx,By,Bz)

% ray_receive=struct('lat',[],'lon',[],'height',[],'hour',[],'minute',[],...
%                    'elevs',[],'bearing',[],'ground_range',[],...
%                    'group_range',[],'absorption',[]);
               
ray_receive = struct('lat',[],'lon',[],'height',[],'hour',[],'minute',[],...
                   'elevs',[],'bearing',[],'ground_range',[],...
                   'group_range',[],'absorption',[], 'electron_density',[],...
                   'phase_path',[], 'refractive_index',[], 'range',[]);

repmat(ray_receive,[1,1]);
origin_lat=raytrace_params.origin_lat;
origin_lon=raytrace_params.origin_lon;
origin_ht=raytrace_params.origin_alt;% ground transmitter location
mode=raytrace_params.mode;
receiver_lat=receiver_location.lat;
receiver_lon=receiver_location.lon;
receiver_alt=receiver_location.alt;
rece_dis=Newton_params.rece_dis;
iteration_num=Newton_params.iteration_num;

% find the midpoint geographic location between the transmitter and
% receiver along great circle path
lat1=origin_lat;
lon1=origin_lon;
lat2=receiver_lat;
lon2=receiver_lon;
lat1_rad = deg2rad(lat1);
lon1_rad = deg2rad(lon1);
lat2_rad = deg2rad(lat2);
lon2_rad = deg2rad(lon2);

% find midpoint latitude
X = cos(lat2_rad) * cos(lon2_rad - lon1_rad);
Y = cos(lat2_rad) * sin(lon2_rad - lon1_rad);
mid_lat_rad = atan2(sin(lat1_rad) + sin(lat2_rad), sqrt((cos(lat1_rad) + X) ^ 2 + Y ^ 2));
% find midpoint longitude
mid_lon_rad = lon1_rad + atan2(Y, cos(lat1_rad) + X);
% convert rad into degree
mid_lat = rad2deg(mid_lat_rad);
mid_lon = rad2deg(mid_lon_rad);

% find quarter latitude
Xquar = cos(mid_lat_rad) * cos(mid_lon_rad - lon1_rad);
Yquar = cos(mid_lat_rad) * sin(mid_lon_rad - lon1_rad);
quarter_lat_rad = atan2(sin(lat1_rad) + sin(mid_lat_rad), sqrt((cos(lat1_rad) + Xquar) ^ 2 + Yquar ^ 2));
% find quarter longitude
quarter_lon_rad = lon1_rad + atan2(Yquar, cos(lat1_rad) + Xquar);
% convert rad into degree
quarter_lat = rad2deg(quarter_lat_rad);
quarter_lon = rad2deg(quarter_lon_rad);

% find quarter latitude
Xeight = cos(quarter_lat_rad) * cos(quarter_lon_rad - lon1_rad);
Yeight = cos(quarter_lat_rad) * sin(quarter_lon_rad - lon1_rad);
eight_lat_rad = atan2(sin(lat1_rad) + sin(quarter_lat_rad), sqrt((cos(lat1_rad) + Xeight) ^ 2 + Yeight ^ 2));
% find quarter longitude
eight_lon_rad = lon1_rad + atan2(Yeight, cos(lat1_rad) + Xeight);
% convert rad into degree
eight_lat = rad2deg(eight_lat_rad);
eight_lon = rad2deg(eight_lon_rad);

wgs84 = wgs84Ellipsoid;

% Determine AER coordinates of first bounce, eg. starting bearing (az) and
% elevation (elev) [Either midpoint or quarter point]
switch mode
    case 1 %1-hop transmission with reflection at 250 km altitude (F-layer)
        H=250e3;
        [az,elev,~] = geodetic2aer(mid_lat,mid_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=1;
    case 2 %1-hop transmission with reflection at 120 km altitude (E-layer)
        H=120e3;
        [az,elev,~] = geodetic2aer(mid_lat,mid_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=1;
    case 3 %2-hop transmission with reflection at 250 km altitude (F-layer)
        H=250e3;
        [az,elev,~] = geodetic2aer(quarter_lat,quarter_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
      nhops=2;
    case 4 %2-hop transmission with reflection at 120 km altitude (E-layer)
        H=120e3;
        [az,elev,~] = geodetic2aer(quarter_lat,quarter_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=2;
    case 5 %1-hop transmission with reflection at 80 km altitude (D-layer)
        H=80e3;
        [az,elev,~] = geodetic2aer(mid_lat,mid_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=1;
    case 6 %2-hop transmission with reflection at 80 km altitude (D-layer)
        H=80e3;
        [az,elev,~] = geodetic2aer(quarter_lat,quarter_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=2;
   
    case 7 %4-hop transmission with reflection at 250 km altitude (F-layer)
        H=250e3;
        [az,elev,~] = geodetic2aer(eight_lat,eight_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=4;
    case 8 %4-hop transmission with reflection at 120 km altitude (E-layer)
        H=120e3;
        [az,elev,~] = geodetic2aer(eight_lat,eight_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=4;
    case 9 %1-hop transmission with reflection at 80 km altitude (D-layer)
        H=80e3;
        [az,elev,~] = geodetic2aer(eight_lat,eight_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=4;
end
fprintf(nhops + "-hop transmission with reflection at " + H/1e3 + " km altitude\n");

if elev < 0
    fprintf("NEGATIVE INITIAL ELEVATION \n")
end
    
initial_elevs=elev;
initial_bear=az;
delta_elevs=raytrace_params.init_delta_elev;
delta_bearing=raytrace_params.init_delta_bear; % initial elevs and bear variation need to adjust with different conditions.

iono_en_grid = Ne;
iono_en_grid_5 = Ne_5;
tic % record time
for Newton_count=1:20  % maximum 20 rounds of Newton's iteration
    if Newton_count>1 % use the smallest distance ray from first round Newton's iteration as new initial
        delta_elevs=rand(1);
        delta_bearing=0.5*rand(1);
        [~,resolve_ray_index]=min(sR(:));
        initial_elevs=elevs_array(resolve_ray_index);
        initial_bear=bearing_array(resolve_ray_index);
    end
    
    elevs_matrix=nan(iteration_num,4);
    bearing_matrix=nan(iteration_num,4);
    sR=nan(iteration_num,4);
    label=nan(1,4);
    
    for m=1:iteration_num
        elevs=[initial_elevs,initial_elevs+delta_elevs];
        bearing_angle=[initial_bear,initial_bear+delta_bearing];% set up elevation and bearing angle
        
        elevs_matrix(m,:)=[initial_elevs,initial_elevs+delta_elevs,initial_elevs,initial_elevs+delta_elevs];
        bearing_matrix(m,:)=[initial_bear,initial_bear,initial_bear+delta_bearing,initial_bear+delta_bearing];
        
        elevs_array=elevs_matrix(:);
        bearing_array=bearing_matrix(:);% record elevation and bearing for updating "inital_elevs" and "initial_bear"
        
        fr=raytrace_params.fr;
        freqs = ones(size(elevs))*fr;
        i=1;
        ray_inf=struct('lat',[],'lon',[],'height',[],'hour',[],'minute',[],'elev',[],'bearing',[],'ground_range',[],...
        'group_range',[],'absorption',[],'electron_density',[]);
   
        repmat(ray_inf,[1,4]);
        lat_receive=nan([1,4]);
        lon_receive=nan([1,4]);
        for bearnumber=1:length(bearing_angle)
            ray_bears =ones(size(elevs))*bearing_angle(bearnumber);
            % number of hops
            tol =  raytrace_params.tol;      % ODE solver tolerance and min max stepsizes
            OX_mode = raytrace_params.OX_mode;
            [ray_data, ray, ~] = ...
                raytrace_3d(origin_lat, origin_lon, origin_ht, elevs, ray_bears, freqs, ...
                OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
                collision_freq, iono_grid_parms, Bx, By, Bz, ...
                geomag_grid_parms);
            assignin('base','ray_data',ray_data);
            assignin('base','ray_path',ray);
            
            for elevsnumber=1:length(elevs)
                if length(ray_data(elevsnumber).ray_label)==1&&ray_data(elevsnumber).ray_label==1
                    start_end = find(ray(elevsnumber).height<=0);
                    endpoint=start_end(2);
                    lat_land=ray(elevsnumber).lat(endpoint);
                    lon_land=ray(elevsnumber).lon(endpoint);
                    absorption_land=ray(elevsnumber).absorption(endpoint);
                elseif length(ray_data(elevsnumber).ray_label)==2&&ray_data(elevsnumber).ray_label(2)==1
                    landpoint=find(ray(elevsnumber).height<=0);
                    endpoint=landpoint(2);
                    absorption_land=ray(elevsnumber).absorption(endpoint);
                else
                    continue
                end
                
                %record selected rays information with cutoff at the potential receiving point.
                ray_inf(i).lat=ray(elevsnumber).lat(1:endpoint);
                ray_inf(i).lon=ray(elevsnumber).lon(1:endpoint);
                ray_inf(i).height=ray(elevsnumber).height(1:endpoint);
                ray_inf(i).elev=ray(elevsnumber).initial_elev;
                ray_inf(i).bearing=ray(elevsnumber).initial_bearing;
                ray_inf(i).absorption=ray(elevsnumber).absorption(1:endpoint);
                ray_inf(i).group_range=ray(elevsnumber).group_range(1:endpoint);
                label(i)=ray_data(elevsnumber).ray_label(end);
                ray_inf(i).electron_density=ray(elevsnumber).electron_density(1:endpoint);
                ray_inf(i).refractive_index=ray(elevsnumber).refractive_index(1:endpoint);
                ray_inf(i).phase_path=ray(elevsnumber).phase_path(1:endpoint);
                lon_receive(i)=ray_inf(i).lon(end);
                lat_receive(i)=ray_inf(i).lat(end);% receiving point latitude and longitude
                i=i+1;
            end
        end
        
        land_index=find(label==1);
        assignin('base','ray_label',label);
        assignin('base','ray_inf',ray_inf);
        
        if length(land_index)<3
            continue
        end
        
        [~,~,newslantRange] = geodetic2aer(receiver_lat,receiver_lon,receiver_alt,lat_receive,lon_receive,receiver_alt,wgs84);
        sR(m,:)=newslantRange;
        [M,receive_ray_index]=min(newslantRange); % find the minimum distance with assumption that ray receiving point at the same altitude as satellite.
        if M>rece_dis % slantRange larger than threshold continue the iteration
            % The iteration equations from James, H. G. (2006) paper
            
            slope1=(lon_receive(1)-lon_receive(3))/(lat_receive(1)-lat_receive(3));
            slope2=(lon_receive(1)-lon_receive(2))/(lat_receive(1)-lat_receive(2));
            delta_elevs=delta_elevs*(receiver_lon-lon_receive(1)-slope1*(receiver_lat-lat_receive(1)))/(lon_receive(2)-lon_receive(1)-slope1*(lat_receive(2)-lat_receive(1)));
            delta_bearing=delta_bearing*(receiver_lon-lon_receive(1)-slope2*(receiver_lat-lat_receive(1)))/(lon_receive(3)-lon_receive(1)-slope2*(lat_receive(3)-lat_receive(1)));
        else
            % slantRange smaller than threshold the received ray is found
            ray_receive.lat=ray_inf(receive_ray_index).lat;% save the received ray information
            ray_receive.lon=ray_inf(receive_ray_index).lon;
            ray_receive.height=ray_inf(receive_ray_index).height;
            ray_receive.elevs=ray_inf(receive_ray_index).elev;
            ray_receive.bearing=ray_inf(receive_ray_index).bearing;
%               ray_receive_V1.loca_num=ray_inf(receive_ray_index).loca_num;
            ray_receive.absorption=ray_inf(receive_ray_index).absorption;
            ray_receive.group_range=ray_inf(receive_ray_index).group_range;
            ray_receive.electron_density=ray_inf(receive_ray_index).electron_density;
            ray_receive.phase_path=ray_inf(receive_ray_index).phase_path;
            ray_receive.refractive_index=ray_inf(receive_ray_index).refractive_index;
            ray_receive.range=M;
%               L=L+1;
            break % slantRange smaller than the threshold break from iteration loop
        end
        assignin('base','m',m);
    end
    
    if length(land_index)<3
            continue
    end
    if M<rece_dis % slantRange smaller than the threshold break from
        break
    end
    
end

% roun_number=Newton_count;

% assignin('base','Newton_count',Newton_count);
% for i=1:length(ray_receive)
%     distance(i)=ray_receive(i).range;
% end
% [~,I]=min(distance);
% assignin('base','receive_ray_index',receive_ray_index);
%  assignin('base','ray_receive',ray_receive);
% assignin('base','distance',distance);
%     loca_num % moniter the process
% end
% select the receiving
% N=1;
% for i=1:length(ray_receive_V1)
% loca_num=ray_receive_V1.loca_num;
% alt_diff=receiver_alt-ray_receive_V1.height(end)*1000;
% if abs(alt_diff)<=rece_dis
%     ray_receive=ray_receive_V1;
%     %         N=N+1;
%     %     end
% else
%     
%     fprintf('\n No ray is received.\n');
% end
% if ~exist('ray_receive','var')
% fprintf('\n No ray is received initial condition need change.\n');
% end

toc
fprintf("\n")
end