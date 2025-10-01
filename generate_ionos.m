% Script meant to be used to generate various ionospheres based on 
% Different R12 values

%%
% Clear console output and memory

clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")
clf

%%
% Define parameters for generation

% ~~~~~~~
% Dummy values to be passed into IONOS constructor - does not affect
% ionosphere generation
date = [0 0 0 0 0]; % Date UTC [YYYY MM DD hh mm]
elevs = [0 1 1];    % Elevation of simulated rays (°)
freq = 10;          % Frequency of simulated rays (MHz)
mode = 1;           % Mode of simulated rays ('X' or 'O' mode)
brk = 0;            % Toggle generation of ray breakdown
% ~~~~~~~

% R12 values to be simulated
R12s = [10, 100];

% Height range, lat range and lon range of ionosphere
% [Start, increment, number of increments]
ht_p  = [0, 4, (24 +1)];   % Height of sim ionosphere (km)
lat_p = [39, 0.5, 4];       % Length of sim ionosphere (° latitude N-S)
lon_p = [-106, 1.0, 10];    % Width of sim ionosphere (° longitude E-W)

% Final toggle to generate ionosphere through IONOS object
% 0 = Do not Generate, 1 = Generate
gen = 0; 

%%
% Core loop - repeat ionosphere generation per R12 value

for sunspot_ind = 1:1:length(R12s)
    R12 = R12s(sunspot_ind);
    fprintf("Simulating R12 value: " + R12 + "\n\n")

    obj = IONS(date, elevs, freq, R12, mode, gen, brk, ht_p, lat_p, lon_p);
end
