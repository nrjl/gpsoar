function [M_frame] = plot_variance(x_grid, y_grid, z_grid, VW_grid, pos_full, box_limits, target, n_isosurf)

figure(2); clf;
h_fig = gcf;
axis equal
set(gca, 'ZDir', 'reverse', 'YDir', 'reverse');

% h_scatter = scatter3(x_grid(:), y_grid(:), z_grid(:), 10, VW_grid(:), 'filled');
if nargin < 6; n_isosurf = 5; end;

minV = min(VW_grid(:));
maxV = max(VW_grid(:));

cmap = colormap;
h_isosurf = zeros(n_isosurf, 1);
V_value = zeros(n_isosurf, 1);

for i = 1:n_isosurf
	V_scale = i/(n_isosurf+1);
	V_value(i) = minV + V_scale*(maxV - minV);
	h_isosurf(i) = patch(isosurface(x_grid, y_grid, z_grid, VW_grid, V_value(i)));
	n_cmap = round(size(cmap, 1)*V_scale);
	set(h_isosurf(i), 'FaceColor', cmap(n_cmap, :), 'EdgeCOlor', 'none',...
		'FaceAlpha', 0.5);
end

camlight
lighting gouraud
hold on;
plot3(pos_full(1,:), pos_full(2,:), pos_full(3,:));
plot3(box_limits(1, [1,2,2,1,1,1,2,2,1,1,1,1,2,2,2,2]), ...
	  box_limits(2, [1,1,1,1,1,2,2,2,2,2,2,1,1,2,2,1]), ...
	  box_limits(3, [1,1,2,2,1,1,1,2,2,1,2,2,2,2,1,1]), ...
	  'Color', [.8, .8, .8], 'LineStyle', '--');
plot3(target(1), target(2), target(3), 'x','Color', [1,0,0], 'MarkerSize', 10);
xlabel('X'); ylabel('Y'); zlabel('Z');
axis tight; view(3);
colorbar('YTickLabel', round(linspace(minV, maxV, 11)*10000)./10000);
if nargout>0
	M_frame = getframe(h_fig, [0, 0, 640, 480]);
end
