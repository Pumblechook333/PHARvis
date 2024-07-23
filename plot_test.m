data = [0,1,100];
x = data;
y = data;
z = data;

k = 1;
kmax = 5;
% plot the rays
while k <= kmax
    
    figure(1)
    pos = get(gcf, 'position');
    pos(3) = pos(3)*1.5;
    pos(4) = pos(4)*1.5;
    set(gcf, 'position', pos)

    set(gca, 'Zlim', [0 400])
    hold on
    
    plot3(x*k, y/k, z*k, 'g')
    
    set(gca, 'XDir','reverse')
    grid on
    hold off
    
    fig = gcf;
    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'uparrow')
      k = k - 1;
      clf
    else
      k = k + 1;
      clf
    end
    
end


