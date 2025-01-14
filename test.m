x = 1:1:10;
y = 1:1:10;

figure(1)
plot(x,y)

UT_str = 'sometime';
filename = "PLOTS/" + UT_str + "2d_raytrace.png";
print(filename, '-dpng')