% Pharlap 3D ground-to-ground honing algorithm

% Coordinates of K2MFF and WWV towers, respectively
k_coord = [40.742018, -74.178975]; w_coord = [40.67583063, -105.038933178];
height = 0; % Assuming from ground
reciever_location = struct('lat', k_coord(1), 'lon', k_coord(2), 'alt', height);
Newton_params = struct('rece_dis', 2594.2027380292675 * 1000, 'iteration_num', 10);

% solver tolerance, minimum step size, maximum step size
tolerance = [1e-7; 0.01; 10];

% Mode to consider:
% 1: 1-hop transmission with reflection at 250 km altitude
% 2: 1-hop transmission with reflection at 120 km altitude
% 3: 2-hop transmission with reflection at 250 km altitude
% 4: 2-hop transmission with reflection at 120 km altitude
mode = 1;   

ide = 0; idb = 1; 
freq = 10;  % Frequency of transmission (10MHz)
OX = 0;     % -1 = Extraordinary, 0 = No field, 1 = Ordinary
nhops = 1;  % Number of hops during transmission

raytrace_params = struct('origin_lat', w_coord(1), 'origin_lon', w_coord(2), 'origin_alt', height,... 
                         'init_delta_elev', ide, 'init_delta_bear', idb, 'fr', freq,...
                         'mode', mode, 'OX_mode', OX, 'nhops', nhops, 'tol', tolerance);
                     
% origin latitude, latitude step, number of latitudes, origin longitude,
% longitude step, number of longitudes, geodetic height, height step,
% number of heights
scale = 100;

max_range = 10000; num_range = 200/scale + 1; range_inc = max_range/(num_range - 1);
height_inc = 3; num_heights = 200/scale;
iono_grid_parms = [w_coord(1), range_inc, num_range, w_coord(2), range_inc,...
                   num_range, height, height_inc, num_heights];
geomag_grid_parms = iono_grid_parms;

% Simulation for 2021/07/01, sunrise
UT = [2021, 7, 1, 10, 32];

% Yearly smoothed monthly median sunspot number
% (Use -1 to search database for date)
R12 = -1;

% generate ionosphere 5 minutes later so that Doppler shift can be calculated
doppler_flag = 1;

% 3d grid of ionospheric electron density, 
% 3d grid of ionospheric electron density 5 minutes later
% 3d grid of effective electron collision frequency (Hz)
% Bx, By Bz: 3d grids of x, y and z components of the geomagnetic field 
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, doppler_flag);

% Placeholders
ray_data = 0; ray_path = 0; ray_label = 0;
ray_inf = 0; m = 0; Newton_count = 0;
receive_ray_index = 0;

ray_recieve = Honing_ground_to_ground(reciever_location, Newton_params,... 
raytrace_params, iono_grid_parms, geomag_grid_parms, iono_pf_grid, ...
iono_pf_grid_5, collision_freq, Bx, By, Bz)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function ray_receive=Honing_ground_to_ground(receiver_location,...
    Newton_params,raytrace_params,iono_grid_parms,geomag_grid_parms,Ne,Ne_5,...
    collision_freq,Bx,By,Bz)

ray_receive=struct('lat',[],'lon',[],'height',[],'hour',[],'minute',[],...
                   'elevs',[],'bearing',[],'ground_range',[],...
                   'group_range',[],'absorption',[]);

repmat(ray_receive,[1,1]);
% speed_of_light = 2.99792458e8;
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

wgs84 = wgs84Ellipsoid;

switch mode
    case 1 %1-hop transmission with reflection at 250 km altitude
        H=200000;
        [az,elev,~] = geodetic2aer(mid_lat,mid_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=1;
    case 2 %1-hop transmission with reflection at 120 km altitude
        H=100000;
        [az,elev,~] = geodetic2aer(mid_lat,mid_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=1;
       
    case 3 %2-hop transmission with reflection at 250 km altitude
        H=200000;
        [az,elev,~] = geodetic2aer(quarter_lat,quarter_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
      nhops=2;
    case 4 %2-hop transmission with reflection at 120 km altitude
        H=100000;
        [az,elev,~] = geodetic2aer(quarter_lat,quarter_lon,H,origin_lat,origin_lon,origin_ht,wgs84);
        nhops=2;
end

% elev_tol=Newton_params.elev_tol;% set up raytracing parameters.
% find initial elevationa and bearing angle
% wgs84 = wgs84Ellipsoid;
% [az,elev,~] = geodetic2aer(receiver_lat,receiver_lon,receiver_alt,origin_lat,origin_lon,origin_ht,wgs84);

% L=1;
% roun_number=nan(1,length(receiver_lat)); % evaluate the iteration efficiency
% for loca_num=1:length(receiver_lat)

initial_elevs=elev;
initial_bear=az;
delta_elevs=raytrace_params.init_delta_elev;
delta_bearing=raytrace_params.init_delta_bear; % initial elevs and bear variation need to adjust with different conditions.

%     if elev(loca_num)<elev_tol %start tratracing procese with initial elevation threshold
%         continue
%     end

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
%             nhops = raytrace_params.nhops; 
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
                    endpoint=find(ray(elevsnumber).height<=0);
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
        if M>rece_dis %slantRange larger than threshold continue the iteration
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

assignin('base','Newton_count',Newton_count);
% for i=1:length(ray_receive)
%     distance(i)=ray_receive(i).range;
% end
% [~,I]=min(distance);
assignin('base','receive_ray_index',receive_ray_index);
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

if ~exist('ray_receive','var')
fprintf('\n No ray is received initial condition need change.\n');
end
toc
end