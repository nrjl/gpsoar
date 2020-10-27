%% Plot flight history field estimate
fig_width = 640; fig_height = 480;
plot_proportion = 1;
plot_plane = false;

if plot_plane
	addpath 'C:\Documents and Settings\n.lawrance\My Documents\Nick\Wing Model\SBXC\'
	addpath 'C:\Documents and Settings\n.lawrance\My Documents\Nick\Wing Model\wing_data\'
	addpath 'C:\Documents and Settings\n.lawrance\My Documents\Nick\Wing Model\graphics\'
	initialise_SBXC
end

% Check type of saved filed
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
set(gca, 'Zdir', 'reverse'); set(gca, 'Ydir', 'reverse'); axis equal;
view(3); hold on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
h_box = plot3(x_limits([1,2,2,1,1,1,2,2,1,1,1,1,2,2,2,2]), ...
    y_limits([1,1,1,1,1,2,2,2,2,2,2,1,1,2,2,1]), ...
    z_limits([1,1,2,2,1,1,1,2,2,1,2,2,2,2,1,1]), ...
    'Color', [.8, .8, .8], 'LineStyle', '--');


% Plot start features
CTrain = squeeze(X_train_full(:,:,1));
WTrain = squeeze(W_train_full(:,:,1));
tTrain = t_train_full(:,1);
cc = max_obs; while CTrain(cc,1)==0; cc=cc-1; end;
CTrain =CTrain(1:cc, :); WTrain =WTrain(1:cc, :); tTrain =tTrain(1:cc, :);

h_train = plot3(CTrain(:,1), CTrain(:,2), CTrain(:,3), 'r+', 'MarkerSize', 3);
h_path = plot3(pos_full(1,1), pos_full(2,1), pos_full(3,1), 'k-');

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
scale = sqrt(max(sum(WW_grid.^2, 2)))/max_wind*0.6;
h_windcone = coneplot(x_grid, y_grid, z_grid, U_grid, V_grid, W_grid, ...
    x_grid, y_grid, z_grid, VW_grid, scale);
set(h_windcone, 'EdgeColor', 'none')
axis tight
h_cbar = colorbar('North');
set(gca, 'xlim', xlim, 'ylim', ylim, 'zlim', zlim);

M3(1) = getframe(h_fig, [0, 0, fig_width, fig_height]);

n_frames = size(t_train_full, 2);
hyper_count2 = 1;

if plot_plane
	X_builder = @(x, att) [zeros(6,1); att(:); x(:)];
	h_plane = SBXC_handle(X_builder(pos_full(:,1), att_full(:,1)), PLANE_AERO);
end

tt = (lookahead*ntf):(lookahead*ntf):tf;
h_time = text(-1000, 500, -25000, 't = 0s');

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
	scale = sqrt(max(sum(WW_grid.^2, 2)))/max_wind*0.6;
	delete(h_windcone);
	h_windcone = coneplot(x_grid, y_grid, z_grid, U_grid, V_grid, W_grid, ...
		x_grid, y_grid, z_grid, VW_grid, scale);
	set(h_windcone, 'EdgeColor', 'none')
% 	axis tight
	
	if plot_plane
		Cr = eye(4); Cr(1:3, 1:3) = calc_Ceb(att_full(:, (ii-1)*replan_points+1));	% Rotation
		Ct = eye(4); Ct(1:3,4) = pos_full(:, (ii-1)*replan_points+1);			% Translation
		set(h_plane, 'Matrix', Ct*Cr)
	end
	
	set(h_time, 'string', sprintf('t = %is', tt(ii)));
	
	M3(ii) = getframe(h_fig, [0, 0, fig_width, fig_height]);
end

%%
if exist('M_title',  'var')
	M3(11:10+ii) = M3(1:ii);
	for jj = 1:10
		M3(jj) = M_title;
	end
end

%% Plot actual field history

if GPt
% Field estimate
h_fig = figure(11); clf; opos = get(h_fig, 'Position');
opos(1:2) = min(opos(1:2), [1920-fig_width, 1100-fig_height]);
set(h_fig, 'Position', [opos(1:2), fig_width, fig_height]);
set(gca, 'Zdir', 'reverse'); set(gca, 'Ydir', 'reverse'); axis equal;
view(3); hold on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
h_box = plot3(x_limits([1,2,2,1,1,1,2,2,1,1,1,1,2,2,2,2]), ...
    y_limits([1,1,1,1,1,2,2,2,2,2,2,1,1,2,2,1]), ...
    z_limits([1,1,2,2,1,1,1,2,2,1,2,2,2,2,1,1]), ...
    'Color', [.8, .8, .8], 'LineStyle', '--');

tt = 0:(lookahead*ntf):tf;

for ii = 1:numel(tt)
	W_true_grid = W_actual([x_grid(:), y_grid(:), z_grid(:)]', ...
		tt(ii)*ones(1, numel(x_grid)));
	scale = sqrt(max(sum(W_true_grid.^2, 1)))/max_wind*0.6;
	h_W_true = coneplot(x_grid, y_grid, z_grid, ...
		reshape(W_true_grid(1,:)', size(x_grid)), ...
		reshape(W_true_grid(2,:)', size(x_grid)), ...
		reshape(W_true_grid(3,:)', size(x_grid)), ...
		x_grid, y_grid, z_grid, scale);
	set(h_W_true, 'EdgeColor', 'none')
% 	h_W_true = quiver3(x_grid, y_grid, z_grid, ...
% 		reshape(W_true_grid(1,:)', size(x_grid)), ...
% 		reshape(W_true_grid(2,:)', size(x_grid)), ...
% 		reshape(W_true_grid(3,:)', size(x_grid)), scale);
	
	M4(ii) = getframe(h_fig, [0, 0, fig_width, fig_height]);
	delete(h_W_true);	
end
end
