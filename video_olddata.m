%% Time-limited estimate video (no future data)
plot_actual = false;		% Plot actual field (need W_actual function)
plot_est = true;			% Plot field estimate
remove_mean = false;			% Remove mean estimate

olddata = true;				% Only plot estimate frum current (old) data
super = false;

%% Super or Ultra version - the choice is yours

if super
	covariance_function = supercov_funs;
	opt_hyper = gfmin;
	n_hyper = supern_hyper;
else
	covariance_function = ultracov_funs;
	opt_hyper = ugfmin;
	n_hyper = ultran_hyper;
end	

lt = exp(opt_hyper(2));
n_nonmeanhyper = sum(n_hyper(1:end-1));
W_steady = zeros(1,3);
W_steady(1:(n_hyper(end)/2)) = opt_hyper(n_nonmeanhyper+(1:2:n_hyper(end)));

% W_steady = [W_mean, 0];

if olddata
	datastop = @(i) i;
else
	datastop = @(i) size(X_train, 1);
end

%%
hfig = figure(3); clf;
haxis = gca;
set(hfig, 'Position', [100, 100, 640, 480]);
daspect([1,1,1]); hold on;

if exist('W_actual', 'var')
	W_grid = W_actual(X_test,t_train(1)*ones(size(X_test,1), 1) );
	max_wind = sqrt(max(sum(W_grid.^2, 2)));
else
	max_wind = max(sqrt(sum(W_train.^2, 2))) - remove_mean*sqrt(sum(W_steady.^2));
end
marker_size = 10;

if plot_actual

if remove_mean
W_grid = (W_grid - ones(size(W_grid, 1), 1)*W_steady);	% REMOVE MEAN WIND AND RESCALE - DODGY
% W_grid = (W_grid - cos_profile3(X_test, cos_layer));	% REMOVE MEAN WIND
end

h_cone_actual = coneplot(x_grid, y_grid, z_grid, ...
	reshape(W_grid(:,1), size(x_grid)), ...
	reshape(W_grid(:,2), size(x_grid)), ...
	reshape(W_grid(:,3), size(x_grid)), ...
	x_grid, y_grid, z_grid);
set(h_cone_actual, 'EdgeColor', 'none', 'FaceColor', [.5 .5 .5]);
end


if plot_est
	[W_est, W_cov] = GPt_predictn_gaussmean(X_train(1:datastop(1),:), t_train(1:datastop(1)), W_train(1:datastop(1),:), ...
		X_test, t_train(1)*ones(size(X_test,1), 1), covariance_function{1}, mean_fun, opt_hyper, n_hyper);
	max_cov = max(W_cov);
	
	if size(W_est, 2) == 2	% Real data
		W_est(:,3) = zeros(size(W_est, 1), 1);
	end
	
	if remove_mean
		W_est = (W_est - ones(size(W_est, 1), 1)*W_steady);	% REMOVE MEAN WIND AND RESCALE - DODGY
		% W_est = (W_est - cos_profile3(X_test, cos_layer));	% REMOVE MEAN WIND
	end
	
	observe_color = square_exp(t_train(1), t_train(1:datastop(i)), log([lt, 1]));
	h_train = scatter3(X_train(1:datastop(i),1), X_train(1:datastop(i),2), X_train(1:datastop(i),3), ...
		marker_size*observe_color, observe_color, 'filled');
	
	est_scale = sqrt(max(sum(W_est.^2, 2)))/max_wind;
	h_cone_est = coneplot(x_grid, y_grid, z_grid, ...
		reshape(W_est(:,1), size(x_grid)), ...
		reshape(W_est(:,2), size(x_grid)), ...
		reshape(W_est(:,3), size(x_grid)), ...
		x_grid, y_grid, z_grid, est_scale, reshape(W_cov, size(x_grid)));
	set(h_cone_est, 'EdgeColor', 'none');
	
else
	h_train = scatter3(X_train(1,1), X_train(1,2), X_train(1,3), ...
		marker_size, 'filled');
end

set(gca, 'ZDir', 'reverse', 'YDir', 'reverse'); axis equal; axis tight;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
M(1) = getframe(hfig);
M = repmat(M(1), 1, n_train);
% fprintf('\n\nFrame 1 of %g', n_train)

for i = 2:n_train
	
	t_test = t_train(i)*ones(size(X_test,1), 1);
	
	if plot_est
	[W_est, W_cov] = GPt_predictn_gaussmean(X_train(1:datastop(i),:), ...
		t_train(1:datastop(i)), W_train(1:datastop(i),:), X_test, ...
		t_test, covariance_function{1}, mean_fun, opt_hyper, n_hyper);
		
	if size(W_est, 2) == 2	% Real data
		W_est(:,3) = zeros(size(W_est, 1), 1);
	end

% 	max_est = sqrt(max(sum(W_grid.^2, 2)));
	
	if remove_mean
		W_est = (W_est - ones(size(W_est, 1), 1)*W_steady);	% REMOVE MEAN WIND
		% W_est = (W_est - cos_profile3(X_test, cos_layer));	% REMOVE MEAN WIND
	end
	
	delete(h_cone_est);
	est_scale = sqrt(max(sum(W_est.^2, 2)))/max_wind;
	h_cone_est = coneplot(x_grid, y_grid, z_grid, ...
		reshape(W_est(:,1), size(x_grid)), ...
		reshape(W_est(:,2), size(x_grid)), ...
		reshape(W_est(:,3), size(x_grid)), ...
		x_grid, y_grid, z_grid, est_scale, reshape(W_cov, size(x_grid)));
	set(h_cone_est, 'EdgeColor', 'none');
	
	observe_color = square_exp(t_train(i), t_train(1:datastop(i)), log([lt, 1]));
	set(h_train, 'XData', X_train(1:datastop(i),1), 'YData', X_train(1:datastop(i),2), ...
		'ZData', X_train(1:datastop(i),3), 'Cdata', observe_color, ...
		'SizeData', marker_size*observe_color);
	else
		set(h_train, 'XData', X_train(1:datastop(i),1), 'YData', X_train(1:datastop(i),2), ...
		'ZData', X_train(1:datastop(i),3));
	end
	
	if plot_actual
	W_grid = W_actual(X_test, t_test);

	if remove_mean
		W_grid = (W_grid - ones(size(W_grid, 1), 1)*W_steady);	% REMOVE MEAN WIND - DODGY
		% W_grid = (W_grid - cos_profile3(X_test, cos_layer));	% REMOVE MEAN WIND
	end
		
	delete(h_cone_actual);
	h_cone_actual = coneplot(x_grid, y_grid, z_grid, ...
		reshape(W_grid(:,1), size(x_grid)), ...
		reshape(W_grid(:,2), size(x_grid)), ...
		reshape(W_grid(:,3), size(x_grid)), ...
		x_grid, y_grid, z_grid);
	set(h_cone_actual, 'EdgeColor', 'none', 'facecolor', [.5 .5 .5]);
	end
	
	
	M(i) = getframe(hfig);
% 	fprintf('\rFrame %g of %g', i, n_train)
end

fprintf('\nAll frames complete; converting to avi...\n');
movie2avi(M, 'GPt_movies\OVERWRITE2.avi', 'fps', 10, 'Compression','none')
fprintf('Finished\n');