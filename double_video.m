%% Plot double video (above below)
fig_width = 640; fig_height = 640;
plot_proportion = 1;

% Check type of saved file
if ~exist('GPt', 'var')
	GPt = false;
end

if exist('save_hypers', 'var')
	retrain = true;
	loghyper2 = save_hypers(1,2:end);
else
	retrain = false;
end

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

xlabel(h_ax(1), 'X'); ylabel(h_ax(1), 'Y'); zlabel(h_ax(1), 'Z');

for ii = 1:2
	h_box(ii) = plot3(h_ax(ii), x_limits([1,2,2,1,1,1,2,2,1,1,1,1,2,2,2,2]), ...
    y_limits([1,1,1,1,1,2,2,2,2,2,2,1,1,2,2,1]), ...
    z_limits([1,1,2,2,1,1,1,2,2,1,2,2,2,2,1,1]), ...
    'Color', [.8, .8, .8], 'LineStyle', '--');
end

%% Plot start features
CTrain = squeeze(X_train_full(:,:,1));
WTrain = squeeze(W_train_full(:,:,1));
tTrain = t_train_full(:,1);
cc = max_obs; while CTrain(cc,1)==0; cc=cc-1; end;
CTrain =CTrain(1:cc, :); WTrain =WTrain(1:cc, :); tTrain =tTrain(1:cc, :);

h_train = plot3(h_ax(1), CTrain(:,1), CTrain(:,2), CTrain(:,3), 'r+', 'MarkerSize', 3);
h_path = plot3(h_ax(1), pos_full(1,1), pos_full(2,1), pos_full(3,1), 'k-');

h_thermals = plot3(h_ax(2), thermals(:,4), thermals(:,5), thermals(:,6), 'o', ...
		'color', [1, 0.6, 0.1], 'markersize', 10, 'LineWidth', 2, ...
		'MarkerFaceColor', [1, 0.6, 0.1]);
set(h_fig, 'currentAxes', h_ax(2));
h_time = text(-100, -150, -450, 't = 0s');

if GPt
	tnow = t0*ones(numel(x_grid), 1);
	[WW_grid, VW_grid] = GPt_predict(CTrain, tTrain, WTrain, ...
		[x_grid(:), y_grid(:), z_grid(:)], tnow, cov_funs{1}, loghyper2);
else
	[WW_grid, VW_grid] = GP_predict(CTrain, WTrain, ...
		[x_grid(:), y_grid(:), z_grid(:)], cov_funs{1}, loghyper2);
end

U_grid = reshape(WW_grid(:,1), size(x_grid));
V_grid = reshape(WW_grid(:,2), size(x_grid));
W_grid = reshape(WW_grid(:,3), size(x_grid));
VW_grid = reshape(VW_grid, size(x_grid));
scale = sqrt(max(sum(WW_grid.^2, 2)))/max_wind*0.5;
h_windcone = coneplot(h_ax(1), x_grid, y_grid, z_grid, U_grid, V_grid, W_grid, ...
    x_grid, y_grid, z_grid, VW_grid, scale);
set(h_windcone, 'EdgeColor', 'none')
axis tight
h_cbar = colorbar('peer', h_ax(1), 'East');
set(h_cbar, 'position', [0.9, 0.05, 0.03, 0.45], 'YAxisLocation', 'right');

tt = 0:(lookahead*ntf):tf;

	W_true_grid = W_actual([x_grid(:), y_grid(:), z_grid(:)]', ...
		tt(ii)*ones(1, numel(x_grid)));
	scale = sqrt(max(sum(W_true_grid.^2, 1)))/max_wind*0.5;
	h_W_true = coneplot(h_ax(2), x_grid, y_grid, z_grid, ...
		reshape(W_true_grid(1,:)', size(x_grid)), ...
		reshape(W_true_grid(2,:)', size(x_grid)), ...
		reshape(W_true_grid(3,:)', size(x_grid)), ...
		x_grid, y_grid, z_grid, scale);
	set(h_W_true, 'EdgeColor', 'none', 'facecolor', 'b')
	set(h_ax, 'xlim', axlim(1,:), 'ylim', axlim(2,:), 'zlim', axlim(3,:));

set(h_ax, 'xlim', axlim(1,:), 'ylim', axlim(2,:), 'zlim', axlim(3,:));
M3(1) = getframe(h_fig, [0, 0, fig_width, fig_height]);

n_frames = size(t_train_full, 2);
hyper_count2 = 1;
tt = (lookahead*ntf):(lookahead*ntf):tf;
%%

for ii = 2:n_frames
	% Plot start features
	CTrain = squeeze(X_train_full(:,:,ii));
	WTrain = squeeze(W_train_full(:,:,ii));
	tTrain = t_train_full(:,ii);
	cc = max_obs; while CTrain(cc,1)==0; cc=cc-1; end;
	CTrain =CTrain(1:cc, :); WTrain =WTrain(1:cc, :); tTrain =tTrain(1:cc, :);
	
	set(h_train, 'XData', CTrain(:,1), 'YData', CTrain(:,2), ...
		'ZData', CTrain(:,3));
	set(h_path, 'XData', pos_full(1,1:(ii-1)*replan_points), ...
		'YData', pos_full(2,1:(ii-1)*replan_points), ...
		'ZData', pos_full(3,1:(ii-1)*replan_points));
	set(h_thermals, 'XData', thermals(:,4)+windspeed(1)*tt(ii), ...
		'YData', thermals(:,5)+windspeed(2)*tt(ii), ...
		'ZData', thermals(:,6)+windspeed(3)*tt(ii));
	set(h_time, 'string', sprintf('t = %is', tt(ii)));
		
	if retrain && (hyper_count2 < hyper_count-2) && (ii == save_hypers(hyper_count2+1,1))
		loghyper2 = save_hypers(hyper_count2+1,2:end);
		hyper_count2 = hyper_count2 + 1;
	end
	
	
	if GPt
		tnow = (ii/n_frames*(tf-t0)+t0)*ones(numel(x_grid), 1);
		[WW_grid, VW_grid] = GPt_predict(CTrain, tTrain, WTrain, ...
			[x_grid(:), y_grid(:), z_grid(:)], tnow, cov_funs{1}, loghyper2);
	else
		[WW_grid, VW_grid] = GP_predict(CTrain, WTrain, ...
			[x_grid(:), y_grid(:), z_grid(:)], cov_funs{1}, loghyper2);
	end
	
	U_grid = reshape(WW_grid(:,1), size(x_grid));
	V_grid = reshape(WW_grid(:,2), size(x_grid));
	W_grid = reshape(WW_grid(:,3), size(x_grid));
	VW_grid = reshape(VW_grid, size(x_grid));
	scale = sqrt(max(sum(WW_grid.^2, 2)))/max_wind*0.5;
	delete(h_windcone);
	h_windcone = coneplot(h_ax(1), x_grid, y_grid, z_grid, U_grid, V_grid, W_grid, ...
		x_grid, y_grid, z_grid, VW_grid, scale);
	set(h_windcone, 'EdgeColor', 'none')
	
	delete(h_W_true);	
	W_true_grid = W_actual([x_grid(:), y_grid(:), z_grid(:)]', ...
		tt(ii)*ones(1, numel(x_grid)));
	scale = sqrt(max(sum(W_true_grid.^2, 1)))/max_wind*0.5;
	h_W_true = coneplot(h_ax(2), x_grid, y_grid, z_grid, ...
		reshape(W_true_grid(1,:)', size(x_grid)), ...
		reshape(W_true_grid(2,:)', size(x_grid)), ...
		reshape(W_true_grid(3,:)', size(x_grid)), ...
		x_grid, y_grid, z_grid, scale);
	set(h_W_true, 'EdgeColor', 'none', 'facecolor', 'b')
	set(h_ax, 'xlim', axlim(1,:), 'ylim', axlim(2,:), 'zlim', axlim(3,:));
	
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
vw = VideoWriter('movies\soaring_DOUBLE.avi', 'Uncompressed AVI');
vw.FrameRate = 10;
open(vw);
writeVideo(vw, M3);
close(vw);
