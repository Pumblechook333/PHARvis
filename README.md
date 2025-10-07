# PHARvis
 Methods for visualizing PHaRLAP raytracing results in MATLAB.

## IONS.m
 Main driver class. Capable of generating ionosphere, performing raytracing, and producing data results for specific simulation parameters.

__Initialization:__
```
% SET necessary parameters for IONS object

date = [2021 7 1 0 0];  % UT Date / Time of simulation [YYYY MM DD hh mm]
elevs = 0:2:90;         % List of ray launch elevations (°)
freq = 10;              % Simulated ray frequency (MHz)
R12 = 57;               % Sunspot Number
mode = 1;               % 1 = O-mode, 0 = X-mode
gen = 0;                % 1 = Generate ionosphere, 0 = Load ionosphere
brk = 0;                % 1 = Perform ray breakdown, 0 = Do not perform

% Call the IONS object
obj = IONS(date, elevs, freq, R12, mode, gen, brk);

% Get the loaded PLASMA FREQUENCY and GEOMAGNETIC grids
iono_series = obj.get_iono_series();

% Get parameters for PLASMA FREQUENCY and GEOMAGNETIC grids 
% used when generating the synthetic ionosphere
iono = obj.get_iono_parms();
iono_grid_parms = iono.iono_grid_parms;
geomag_grid_parms = iono.geomag_grid_parms;
iono_height = iono.iono_height;

% Get launch elevation, frequency, and datetime of simulated ray
g = obj.get_gen_params();
elevs = g.elevs;
freqs = g.freqs;
UT = g.UT;

% Get coordinates and height above sea level of simulated 
% ray origin station and reciever station
c = obj.get_coords();
origin_lat = c.origin_lat;
origin_lon = c.origin_lon;
origin_ht = c.origin_ht;
RX_coord = c.RX_coord;
TX_coord = c.TX_coord;

% Get launch bearing of simulated rays
ray_bears = obj.get_bearing();
```

## Example figures:

![image](https://github.com/user-attachments/assets/1ba732ba-6a95-4ca9-b509-5b3559bb1d50)

![image](https://github.com/user-attachments/assets/f74e276b-35e5-4c06-a94a-526b52f73b77)

![image](https://github.com/user-attachments/assets/860842b1-ca98-4f0b-bc81-07f8e993dcb3)


