clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n")

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% GET necessary vars

date = [2021 7 1 0 0];

el_start = 0;
% el_inc = 0.2;
% el_stop = 50;
el_inc = 5;
el_stop = 90;
elevs = el_start : el_inc : el_stop;

freq = 10;
R12_sel = [-1, 25, 50, 100, 200];
R12 = R12_sel(5);
gen = 1; % 0 = no gen, 1 = gen

elevs_string = " || Initial Elevations: " ...
               + el_start + ":" + el_inc + ":" + el_stop;
r12_string = " || R12: " + R12;

obj = IONS(date, elevs, freq, R12, gen);

props = ["phase_path", "geometric_distance"];
obj.ray_props(props)

