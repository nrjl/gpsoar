function [h0 h1] = lighten_axes(h0)

% Lighten gridlines
% h0 = gca;					% Original (set to lighter color)
% xlim = get(h0, 'xlim'); ylim = get(h0, 'ylim'); zlim = get(h0, 'zlim'); 
h1 = copyobj(h0,get(h0, 'parent'));		% New (black but no gridlines)

set(get(h0, 'children'), 'visible', 'off');

linkaxes([h1, h0]);

set(h1, 'color', 'none', 'xgrid', 'off', 'ygrid','off');
set(h0, 'Xcolor', [.4 .4 .4]);
set(h0, 'Ycolor', [.4 .4 .4]);

set(get(h0, 'parent'), 'currentaxes', h0);

