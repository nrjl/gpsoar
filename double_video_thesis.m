%% Plot double video (above below)
% plot_result;
fig_width = 640; fig_height = 640;
plot_proportion = 1;

retrain = false;
loghyper2 = loghyper; %save_hypers(1,2:end);

% Field estimate
h_fig = figure(10); clf; opos = get(h_fig, 'Position');
opos(1:2) = min(opos(1:2), [1920-fig_width, 1100-fig_height]);
set(h_fig, 'Position', [opos(1:2), fig_width, fig_height]);

h_ax(1) = axes('position', [0.1, 0.0, 0.8, 0.6]);
h_ax(2) = axes('position', [0.1, 0.45, 0.8, 0.6]);

set(h_ax, 'Zdir', 'reverse'); set(h_ax, 'Ydir', 'reverse');
set(h_ax, 'DataAspectRatio', [1 1 1]);
set(h_ax, 'NextPlot', 'add');
view(h_ax(1), 3); view(h_ax(2), 3);

xlabel(h_ax(1), '\it x_i \rm (m)'); ylabel(h_ax(1), '\it y_i \rm (m)');
zlabel(h_ax(1), '\it z_i \rm (m)');

for ii = 1:2
	h_box(ii) = plot3(h_ax(ii), x_limits([1,2,2,1,1,1,2,2,1,1,1,1,2,2,2,2]), ...
    y_limits([1,1,1,1,1,2,2,2,2,2,2,1,1,2,2,1]), ...
    z_limits([1,1,2,2,1,1,1,2,2,1,2,2,2,2,1,1]), ...
    'Color', [.8, .8, .8], 'LineStyle', '--');
end

set(h_ax, 'xlim', axlim(1,:), 'ylim', axlim(2,:), 'zlim', axlim(3,:));

%% Plot start features
CTrain = squeeze(X_train_full(:,:,1));
WTrain = squeeze(W_train_full(:,:,1));
tTrain = t_train_full(:,1);
cc = max_obs; while CTrain(cc,1)==0; cc=cc-1; end;
CTrain =CTrain(1:cc, :); WTrain =WTrain(1:cc, :); tTrain =tTrain(1:cc);

h_train = plot3(h_ax(1), CTrain(:,1), CTrain(:,2), CTrain(:,3), 'ko', ...
	'MarkerSize', 3, 'markerfacecolor', 'k');
h_path = plot3(h_ax(1), pos_full(1,1), pos_full(2,1), pos_full(3,1),...
	'-', 'linewidth', 0.5, 'color', .4*[1 1 1]);
h_newpath = plot3(h_ax(1), pos_full(1,1), pos_full(2,1), pos_full(3,1),...
	'k-', 'linewidth', 1.5);

h_thermals =zeros(size(thermals, 1), 2);
[xcyl, ycyl, zcyl] = cylinder(0.75); zcyl = zcyl-.5;

for ii = 1:size(thermals, 1)
	h_thermals(ii, 1) = plot3(h_ax(2), thermals(ii,4), thermals(ii,5), ...
		thermals(ii,6), 'o', 'color', [1, 0.6, 0.1], 'markersize', 10, ...
		'MarkerFaceColor', [1, 0.6, 0.1]);
	h_thermals(ii,2) = surf(xcyl*thermals(ii,2)+thermals(ii,4), ...
		ycyl*thermals(ii,2)+thermals(ii,5), ...
		zcyl*thermals(ii,2)+thermals(ii,6), 'facecolor', [1, 0.6, 0.1], ...
		'edgecolor', 'none', 'facealpha', 0.2);
end
	
set(h_fig, 'currentAxes', h_ax(2));
h_time = text(-100, -150, -450, 't = 0s');

tnow = t0*ones(numel(x_grid), 1);
% [WW_grid, VW_grid] = GP_predict(CTrain, WTrain, ...
% 	[x_grid(:), y_grid(:), z_grid(:)], tnow, cov_funs{1}, loghyper2);
[WW_grid, VW_grid] = GPt_predict(CTrain, tTrain, WTrain, ...
	[x_grid(:), y_grid(:), z_grid(:)], tnow, cov_funs{1}, loghyper2);
	

U_grid = reshape(WW_grid(:,1), size(x_grid));
V_grid = reshape(WW_grid(:,2), size(x_grid));
W_grid = reshape(WW_grid(:,3), size(x_grid));
VW_grid = reshape(VW_grid, size(x_grid));
scale = sqrt(max(sum(WW_grid.^2, 2)))/max_wind*0.5;
h_windcone = coneplot(h_ax(1), x_grid, y_grid, z_grid, U_grid, V_grid, W_grid, ...
    x_grid, y_grid, z_grid, VW_grid, scale);
set(h_windcone, 'EdgeColor', 'none', 'facealpha', 0.5);

h_cbar = colorbar('peer', h_ax(1), 'East');
set(h_cbar, 'position', [0.9, 0.05, 0.03, 0.45], 'YAxisLocation', 'right');

tt = 0:(lookahead*ntf):tf;

W_true_grid = W_actual([x_grid(:), y_grid(:), z_grid(:)]', ...
t0);
scale = sqrt(max(sum(W_true_grid.^2, 1)))/max_wind*0.5;
h_W_true = coneplot(h_ax(2), x_grid, y_grid, z_grid, ...
	reshape(W_true_grid(1,:)', size(x_grid)), ...
	reshape(W_true_grid(2,:)', size(x_grid)), ...
	reshape(W_true_grid(3,:)', size(x_grid)), ...
	x_grid, y_grid, z_grid, scale);
set(h_W_true, 'EdgeColor', 'none', 'facecolor', 'b')

h_slice = slice(h_ax(1), x_grid, y_grid, z_grid, VW_grid, max(x_grid(:)), ...
	min(y_grid(:)), max(z_grid(:)));
set(h_slice, 'edgecolor', 'none', 'facecolor', 'interp', 'facealpha', 0.2)

set(h_ax, 'xlim', axlim(1,:), 'ylim', axlim(2,:), 'zlim', axlim(3,:));
M3(1) = getframe(h_fig, [0, 0, fig_width, fig_height]);


%% Time info
actual_fps = 0.8; % Frequency of plotted frames

n_frames = (tf-t0)*actual_fps +1;
tvec = t0:(1/actual_fps):tf;

hyper_count2 = 1;
tt = (lookahead*ntf):(lookahead*ntf):tf;
%%

for ii = 2:n_frames
	tnow = ii/actual_fps;
	
	if ~rem(tnow, replan)
		% Get new feature set and training data
		CTrain = squeeze(X_train_full(:,:,tnow/replan));
		WTrain = squeeze(W_train_full(:,:,tnow/replan));
		tTrain = t_train_full(:,tnow/replan);
		cc = max_obs; while CTrain(cc,1)==0; cc=cc-1; end;
		CTrain =CTrain(1:cc, :); WTrain =WTrain(1:cc, :); tTrain =tTrain(1:cc, :);
		
		set(h_train, 'XData', CTrain(:,1), 'YData', CTrain(:,2), ...
		'ZData', CTrain(:,3));
	end
	
	set(h_path, 'XData', pos_full(1,1:tnow/dt), ...
		'YData', pos_full(2,1:tnow/dt), ...
		'ZData', pos_full(3,1:tnow/dt));
	iimin = max(tnow/dt-20/dt, 1);
	set(h_newpath, 'XData', pos_full(1,iimin:tnow/dt), ...
		'YData', pos_full(2,iimin:tnow/dt), ...
		'ZData', pos_full(3,iimin:tnow/dt));
	
	set(h_thermals(:,1), 'XData', thermals(:,4)+windspeed(1)*tnow, ...
		'YData', thermals(:,5)+windspeed(2)*tnow, ...
		'ZData', thermals(:,6)+windspeed(3)*tnow);
	for jj = 1:size(thermals, 1)
		set(h_thermals(jj,2), 'XData', xcyl*thermals(jj,2)+thermals(jj,4)+windspeed(1)*tnow, ...
		'YData', ycyl*thermals(jj,2)+thermals(jj,5)+windspeed(2)*tnow, ...
		'ZData', zcyl*thermals(jj,2)+thermals(jj,6)+windspeed(3)*tnow);
	end
	
	set(h_time, 'string', sprintf('t = %0.1fs', tnow));
		
		
	nowhyper = sum(save_hypers(:,1) < tnow/replan);
	
	loghyper = save_hypers(nowhyper, 2:end);
	
	tnowvec = tnow*ones(numel(x_grid), 1);
	[WW_grid, VW_grid] = GPt_predict(CTrain, tTrain, WTrain, ...
		[x_grid(:), y_grid(:), z_grid(:)], tnowvec, cov_funs{1}, loghyper2);
	
	U_grid = reshape(WW_grid(:,1), size(x_grid));
	V_grid = reshape(WW_grid(:,2), size(x_grid));
	W_grid = reshape(WW_grid(:,3), size(x_grid));
	VW_grid = reshape(VW_grid, size(x_grid));
	scale = sqrt(max(sum(WW_grid.^2, 2)))/max_wind*0.5;
	delete(h_windcone);
	h_windcone = coneplot(h_ax(1), x_grid, y_grid, z_grid, U_grid, V_grid, W_grid, ...
		x_grid, y_grid, z_grid, VW_grid, scale);
	set(h_windcone, 'EdgeColor', 'none', 'facealpha', 0.5);
	
	delete(h_W_true);	
	W_true_grid = W_actual([x_grid(:), y_grid(:), z_grid(:)]', ...
		tnow);
	scale = sqrt(max(sum(W_true_grid.^2, 1)))/max_wind*0.5;
	h_W_true = coneplot(h_ax(2), x_grid, y_grid, z_grid, ...
		reshape(W_true_grid(1,:)', size(x_grid)), ...
		reshape(W_true_grid(2,:)', size(x_grid)), ...
		reshape(W_true_grid(3,:)', size(x_grid)), ...
		x_grid, y_grid, z_grid, scale);
	set(h_W_true, 'EdgeColor', 'none', 'facecolor', 'b')
	set(h_ax, 'xlim', axlim(1,:), 'ylim', axlim(2,:), 'zlim', axlim(3,:));
	
	
	delete(h_slice);
	h_slice = slice(h_ax(1), x_grid, y_grid, z_grid, VW_grid, max(x_grid(:)), ...
		min(y_grid(:)), max(z_grid(:)));
	set(h_slice, 'edgecolor', 'none', 'facecolor', 'interp', 'facealpha', 0.2)
	
	M3(ii) = getframe(h_fig, [0, 0, fig_width, fig_height]);
	
	
end

%%
if exist('M_title',  'var')
	M3(11:10+ii) = M3(1:ii);
	for jj = 1:10
		M3(jj) = M_title;
	end
end

%%
% movie2avi(M3, 'movies\soaring_DOUBLE2.avi', 'fps', 15, 'Compression','none')
