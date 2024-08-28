classdef IONS
    properties(Constant)
        % 25, 50, 100, 200
        % R12 = -1;            % Autoselect sunspot number
        doppler_flag = 0;    % Not interested in Doppler shift
        TX_coord = [40.67583063, -105.038933178];   % WWV Station
        RX_coord = [40.742018, -74.178975];         % K2MFF Station
        
        % Bottom left corner of grid (lat, lon) in [deg]
        grid_corner = [39, -106];
%         grid_corner = [40, -(105+(48/60))];
%         grid_corner = [36, -110];

    end
    
    properties
        % Gen Params
        date = [2021 7 1 0 0];
        UT
        elevs
        freqs           % Frequency (MHz)
        R12
        mode
        nhops_max = 4;
        
        % Coords
        origin_lat      % latitude of the start point of rays
        origin_lon      % longitude of the start point of rays
        origin_ht       % altitude of the start point of rays
        receiver_lat 
        receiver_lon
        
        % Bearing
        ray_bears
        
        % Iono Parms
        iono_grid_parms
        geomag_grid_parms
        iono_height
        
        % gl_iono
        iono_series = struct();
        
        % Ray Breakdown
        per_tot = 0;
        per_hop = struct();
        
    end
    
    methods(Static)      
        function in_range = chk_dist(range, ray_N, wgs84)
            
            if ~exist('wgs84', 'var')
                wgs84 = wgs84Ellipsoid('km');
            end
            
            hrange = range;
            
            if ray_N.height(end) < hrange
                fLat = ray_N.lat(end);
                fLon = ray_N.lon(end);
                RXLat = IONS.RX_coord(1);
                RXLon = IONS.RX_coord(2);

                d = distance(fLat, fLon, RXLat, RXLon, wgs84);

                % Check if endpoint of ray is within ground range
                if d < range
                    in_range = 1;
                else
                    in_range = 0;
                end

            else
                in_range = -1;
            end
            
        end
        
    end
    
    methods
        % Initialization
        function self = IONS(date, elevs, freq, R12, mode, gen, brk)
            arguments
                date = [2021 7 1 0 0]
                elevs = 0:5:90
                freq = 10
                R12 = -1
                mode = 1
                gen = 0
                brk = false
            end
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % General Parameters
            
            self.R12 = R12;
            self.mode = mode;
            
            if mode == 1
                disp_mode = 'O';
            else
                disp_mode = 'X';
            end 
            fprintf(disp_mode + "-mode, " + R12 + " R12 \n");

            % Generate an Ionosphere for each hour
            self.date = date;           % yyyy m d 0 0
            UT = repmat(self.date, [24,1]);
            for i=1:1:24
                UT(i,4) = i-1;
            end
            self.UT = UT;

            self.elevs = elevs;
            self.freqs = ones(size(self.elevs))*freq;   % frequency (MHz)
            
            self = self.coords();
            self = self.bearing();
            self = self.iono_parms();
            self = self.gl_iono(gen);
            
            if brk 
                tmp_props = [["ground_range", "ray_data"]];
                tmp_rps = self.ray_props(tmp_props);
                self = self.ray_breakdown(tmp_props, tmp_rps);
            end
        end
        
        % Setters
        function self = coords(self)
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Coordinates of Towers

            self.origin_lat = self.TX_coord(1);             % latitude of the start point of rays
            self.origin_lon = self.TX_coord(2);            % longitude of the start point of rays
            self.origin_ht = 0.0;                % altitude of the start point of rays
            self.receiver_lat = self.RX_coord(1); 
            self.receiver_lon = self.RX_coord(2);

        end
        
        function self = bearing(self) 
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Determine Bearing

            % find the midpoint geographic location between the transmitter and
            % receiver along great circle path
            lat1=self.origin_lat;
            lon1=self.origin_lon;
            lat2=self.receiver_lat;
            lon2=self.receiver_lon;
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

            wgs84 = wgs84Ellipsoid;

            H = 250e3; % (F-layer)
            [az,~,~] = geodetic2aer(mid_lat,mid_lon,H,self.origin_lat,self.origin_lon,self.origin_ht,wgs84);
            aim = 0;
            az = az + aim;
            self.ray_bears = ones(size(self.elevs))*az ; % initial bearing of rays            
            fprintf("Bearing of: " + az + "° \n");
            
        end
        
        function self = iono_parms(self, ht_p, lat_p, lon_p)
            arguments
                self
                
                ht_p  = [self.origin_ht, 4, (125 +1)];
                lat_p = [self.grid_corner(1), 0.5, 4];
                lon_p = [self.grid_corner(2), 1.0, 36];
            end
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Ionospheric Grid Parameters

            ht_start = ht_p(1);          % start height for ionospheric grid (km)
            ht_inc = ht_p(2);             % height increment (km)
            num_ht = ht_p(3);     
            iono_ht = ht_inc * (num_ht - 1);
            self.iono_height = iono_ht;
            fprintf("Ceiling at: " + iono_ht + "km \n");
            

            lat_start = lat_p(1);
            lat_inc = lat_p(2);
            lines_lat = lat_p(3);
            num_lat = lines_lat/lat_inc +1;

            lon_start= lon_p(1);
            lon_inc = lon_p(2);
            lines_lon = lon_p(3);
            num_lon = lines_lon/lon_inc +1;

            iono_grid_p = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
                  ht_start, ht_inc, num_ht];
            self.iono_grid_parms = iono_grid_p;

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

            geomag_grid_p = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
                  B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];
            self.geomag_grid_parms = geomag_grid_p;
        end
        
        function self = gl_iono(self, gen, rsta, rinc, rsto)
            arguments
                self
                
                gen = 0;
                rsta = 1;
                rinc = 1;
                rsto = 24;
            end
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Generate / Load Ionosphere
            
            append = "_" + self.iono_height + "km_series_" + self.R12;
            name = 'ionosphere' + append;

            dir = name + '/';
            fname = name + '.mat';

            PARMSPATH = dir + 'params_' + fname;

            if gen == 1
                fprintf('Generating ionospheric and geomag grids... \n')
                if not(isfolder(dir))
                    mkdir(dir)
                end
                
                iono_grid_parms = self.iono_grid_parms; %#ok<*PROPLC>
                geomag_grid_parms = self.geomag_grid_parms;
                UT = self.UT;
                
                % Save the ionospheric grid parameters for later reference
                save(PARMSPATH, 'iono_grid_parms', 'geomag_grid_parms', 'UT')
                for i=rsta:rinc:rsto
                    fprintf('Hour:' + string(i))
                    fprintf('\n')

                    tic

                    [iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
                    gen_iono_grid_3d(self.UT(i,:), self.R12, self.iono_grid_parms, self.geomag_grid_parms, self.doppler_flag);

                    IONOPATH = dir + i + "_" + fname;
                    save(IONOPATH, 'iono_pf_grid', 'iono_pf_grid_5' ,...
                        'collision_freq', 'Bx', 'By', 'Bz')

                    data = struct('iono_pf_grid', iono_pf_grid, 'iono_pf_grid_5',...
                    iono_pf_grid_5, 'collision_freq', collision_freq,...
                    'Bx', Bx, 'By', By, 'Bz', Bz);        

                    field = 'i' + string(i);
                    self.iono_series.(field) = data;

                    toc
                end     
            else
                % Load iono/geomag parameters
                load = matfile(PARMSPATH);
                self.iono_grid_parms = load.iono_grid_parms;
                self.geomag_grid_parms = load.geomag_grid_parms;
                self.UT = load.UT;
                self.iono_height = self.iono_grid_parms(8) * (self.iono_grid_parms(9) - 1);
                
                if ~(self.iono_height == self.iono_height)
                    print("WARNING: IONS.m iono height does not match loaded iono height. \n")
                end
                
                fprintf('Loading ionospheric and geomag grids... \n')
                for i=rsta:rinc:rsto       
                    IONOPATH = dir + i + "_" + fname;
                    load = matfile(IONOPATH);

                    iono_pf_grid = load.iono_pf_grid;
                    iono_pf_grid_5 = load.iono_pf_grid_5;
                    collision_freq = load.collision_freq;
                    Bx = load.Bx;
                    By = load.By;
                    Bz = load.Bz;

                    data = struct('iono_pf_grid', iono_pf_grid, 'iono_pf_grid_5',...
                    iono_pf_grid_5, 'collision_freq', collision_freq,...
                    'Bx', Bx, 'By', By, 'Bz', Bz);        

                    field = 'i' + string(i);
                    self.iono_series.(field) = data;
                end
            end

            fprintf('\n')
            
        end
        
        function res = ray_props(self, props)
            arguments
                self
                
                props
            end
            
            sz = size(props);
            nprops = sz(1);
            
            rsta = 1;
            rinc = 1;
            rsto = 24;

            tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
            num_elevs = length(self.elevs); 
            
            return_properties = struct();
            for p = 1:nprops
                return_properties.(props(p)) = [];
            end

            range = 100;
            wgs84 = wgs84Ellipsoid('km');
            
            for nhops = 1:1:self.nhops_max
                hop_field = "hop_" + nhops;
                
                tmp = zeros(rsto,num_elevs);
                arrs = repmat(tmp, nprops, 1);
                
                for i = rsta:rinc:rsto
                    hour_field = 'i' + string(i);

                    iono_pf_grid = self.iono_series.(hour_field).iono_pf_grid;
                    iono_pf_grid_5 = self.iono_series.(hour_field).iono_pf_grid_5;
                    collision_freq = self.iono_series.(hour_field).collision_freq;
                    Bx = self.iono_series.(hour_field).Bx;
                    By = self.iono_series.(hour_field).By;
                    Bz = self.iono_series.(hour_field).Bz;

                    % convert plasma frequency grid to  electron density in electrons/cm^3
                    iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
                    iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

                    [ray_data_N, ray_N, ~] = ...
                      raytrace_3d(self.origin_lat, self.origin_lon, self.origin_ht, self.elevs, self.ray_bears, self.freqs, ...
                                  self.mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
                                  collision_freq, self.iono_grid_parms, Bx, By, Bz, ...
                                  self.geomag_grid_parms);
                    
                    tmp = zeros(1,num_elevs);
                    ray_props = repmat(tmp, nprops, 1);
                    
                    for elev = 1:1:num_elevs
                        grounded_ray = ray_data_N(elev).ray_label == 1;
                        if grounded_ray
                            in_range = self.chk_dist(range, ray_N(elev), wgs84) == 1;
                            if in_range
                                for p = 1:nprops
                                    if props(p,2) == "ray"
                                        ray_props(p, elev) = ...
                                                ray_N(elev).(props(p))(end);
                                    elseif props(p,2) == "ray_max"
                                        ray_props(p, elev) = ...
                                                max(ray_N(elev).(props(p)));
                                    else
                                        ray_props(p, elev) = ...
                                                  sum(ray_data_N(elev).(props(p)));
                                    end
                                end
                            end
                        end
                    end
                    
                    for p = 1:nprops
                        arrs(i+((p-1)*rsto),:) = ray_props(p, :);
                    end

                end
                
                for p = 1:nprops
                    return_properties.(props(p)).(hop_field) = ...
                    arrs((1 + (rsto*(p-1))):(rsto + rsto*(p-1)), :);
                end

                % finished ray tracing with this ionosphere so clear it out of memory
                clear raytrace_3d
            end
            
            res = return_properties;
            
        end
        
        function self = ray_breakdown(self, properties, return_properties)
            prop1 = properties(1,1);
            samp_arr = return_properties.(prop1);
            num_rays = length(self.elevs);
            
            breakdown = struct();
            
            rsto = 24;
            per_hour_tot = zeros(rsto,1);
            for nhops = 1:1:self.nhops_max
                per_hour = zeros(rsto,1);
                hop_field = "hop_" + nhops;
                
                for hour = 1:1:rsto
                    hour_slice = samp_arr.(hop_field)(hour, :);
                    percent = nnz(hour_slice) / num_rays;
                    per_hour(hour) = percent;
                end
                
                per_hour_tot = per_hour_tot + per_hour;
                breakdown.(hop_field) = per_hour;
                
            end
            self.per_tot = per_hour_tot;
            self.per_hop = breakdown;
            
        end
        
        % Getters
        function res = get_gen_params(self)
            res = struct('date', self.date, 'UT', self.UT, 'elevs', self.elevs, ...
                  'freqs', self.freqs, 'R12', self.R12, 'doppler_flag', self.doppler_flag);
        end
        
        function res = get_coords(self)
            res = struct('TX_coord', self.TX_coord, 'RX_coord', self.RX_coord, 'origin_lat', ...
            self.origin_lat, 'origin_lon', self.origin_lon, 'origin_ht', self.origin_ht, ...
            'receiver_lat', self.receiver_lat, 'receiver_lon', self.receiver_lon);
        end
        
        function res = get_bearing(self)
            res = self.ray_bears;
        end
        
        function res = get_iono_parms(self)
            res = struct('iono_grid_parms', self.iono_grid_parms, ...
                'geomag_grid_parms', self.geomag_grid_parms, ...
                'iono_height', self.iono_height);
        end
        
        function res = get_iono_series(self)
            res = self.iono_series;
        end
        
        function res = get_ray_breakdown(self)
            res = struct("per_tot", self.per_tot, "per_hop", self.per_hop);
        end
    end
end