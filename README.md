# PHARvis
 Methods for visualizing pharlap raytracing results in MATLAB

## IONS.m
 Class containing all necessary elements of Ionosphere

__Initialization:__
```
date = [2021 7 1 0 0];

el_start = 0;
el_inc = 0.2;
el_stop = 50;
elevs = el_start : el_inc : el_stop;

freq = 10;
R12 = 200;

obj = IONS(date, elevs, freq, R12, gen);
```

__Get ray properties:__
```
props = [["ground_range", "ray_data"];
         ["geometric_path_length", "ray_data"]];
rps = obj.ray_props(props);

ground_range = rps.(props(1));
geometric_path_length = rps.(props(2));
```

## Figures:

![image](https://github.com/user-attachments/assets/1ba732ba-6a95-4ca9-b509-5b3559bb1d50)

![image](https://github.com/user-attachments/assets/f74e276b-35e5-4c06-a94a-526b52f73b77)

![image](https://github.com/user-attachments/assets/860842b1-ca98-4f0b-bc81-07f8e993dcb3)


