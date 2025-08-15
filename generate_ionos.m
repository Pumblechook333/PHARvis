% Script meant to be used to generate various ionospheres based on 
% Different R12 values

clear
clc
fprintf("~~~~~ " + mfilename + " ~~~~~ \n\n")

clf

% 2023 10 10-20
% year = 2023;
% month = 10;
% days = 10:1:20;
% R12s = [132, 145, 137, 113, 111, 104, 97, 74, 66, 57, 67]; % 2023 10 10-20

% 2024 04 05-15
year = 2024;
month = 04;
days = 5:1:15;
R12s = [78, 80, 84, 77, 55, 55, 90, 86, 120, 148, 171];


% for day = 1:1:11
for day = 5:1:5
    % for r12_i = 1:1:r12_max
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % GET necessary vars

    date = [year month days(day) 0 0];

    fprintf(mat2str(date))

    el_start = 0;
    el_inc = 0.1;
    el_stop = 50;
    elevs = el_start : el_inc : el_stop;

    freq = 10;
    R12 = R12s(day);
    mode = 'O';
    gen = 1; % 0 = no gen, 1 = gen
    
    obj = IONS(date, elevs, freq, R12, mode, gen);

end
