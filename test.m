% rng('default')
% fig = figure();
% tcl = tiledlayout(fig,2,2);
% n = 4;  % number of heatmaps
% h = gobjects(n,1); 
% for i = 1:n
%     ax = nexttile(tcl); 
%     h(i) = heatmap(rand(5)*randi(5),'ColorbarVisible','off');
% end
% % Equate color limits in all heatmaps
% colorLims = vertcat(h.ColorLimits);
% globalColorLim = [min(colorLims(:,1)), max(colorLims(:,2))];
% set(h, 'ColorLimits', globalColorLim)
% % Create global colorbar that uses the global color limits
% ax = axes(tcl,'visible','off','Colormap',h(1).Colormap,'CLim',globalColorLim);
% cb = colorbar(ax);
% set(cb, 'Position', [0.85 0.15 0.05 0.7]);

rng('default')
fig = figure();
tcl = tiledlayout(fig,2,2);
n = 4;  % number of heatmaps
h = gobjects(n,1); 
for i = 1:n
    ax = nexttile(tcl); 
    h(i) = heatmap(rand(5)*randi(5),'ColorbarVisible','off');
end

% Equate color limits in all heatmaps
colorLims = vertcat(h.ColorLimits);
globalColorLim = [min(colorLims(:,1)), max(colorLims(:,2))];
set(h, 'ColorLimits', globalColorLim)

% Create a new axes for the colorbar, ensuring it's visible
ax_cb = axes('Position',[0.95 0.15 0.05 0.7]);
colormap(ax_cb, h(1).Colormap);
colorbar(ax_cb);
% caxis(ax_cb, globalColorLim);