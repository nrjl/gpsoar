%% Setup for running GP over flight wind data
clear variables
close all

plot_on = 1;
%% GENERATE OR RETRIEVE FULL SET OF TRAINING DATA

thermal_props = [4, 100, 3, 0, 0, -200];
W_steady = [1, 0, -0.5];
cos_layer = [-200, -100, 3, 0, 1, 3];  %[zl, zu, vl, vu, dW, ds]
cos_layer2 = [-300, -250, -3, 0, 1, 3];  %[zl, zu, vl, vu, dW, ds]

% W_actual = @(X, t) thermal_field((X - t*W_steady), thermal_props);

% t_shift = @(XX, tt, pp) (XX - tt*pp);

W_actual = @(X, t) multi_field3(X, ...
	(@(XX, p) thermal_field((XX-(t*W_steady)'), p)), thermal_props);%, ...
% 	@cos_profile3, cos_layer, ...
% 	@cos_profile3, cos_layer2);
% 	@(XX, p) linear_wind3(XX, p, zeros(3)), [W_steady(1:2)'; 0]);
% 	@(XX, pp) shift_transform(XX, t, pp, @(XX, tt, pp) (XX - tt*pp), W_steady, @thermal_field), thermal_props, ...
% W_actual = @(X, t) cos_profile3(X, cos_layer);



sigma_s = 0.05;


%% TRAINING DATA
% These are the training data points, i.e. the points that you have sampled
% and from which you wish to reconstruct the field.

% -- Helical path -- %
ri = 50;
ro = 100;
zb = -100; zt = -300;
n_cycles = 3;
n_samples = 200;
mean_speed = 15;

% Single helix (climbing)
% [x_train, y_train, z_train] = ...
%  	helical_path(ri, ro, zb, zt, n_cycles, n_samples);

% Double helix (climbing then sinking)
[x_train1, y_train1, z_train1] = ...
 	helical_path(ri, ro, zb, zt, n_cycles, floor((n_samples+1)/2));
[x_train2, y_train2, z_train2] = ...
 	helical_path(ri, ro, zt, zb, n_cycles, ceil((n_samples+1)/2));
x_train = [x_train1(:); x_train2(2:end)'];
y_train = [y_train1(:); y_train2(2:end)']; 
z_train = [z_train1(:); z_train2(2:end)']; 

X_train = [x_train(:), y_train(:), z_train(:)];
path_length = sqrt(sum((X_train(2:end,:) - X_train(1:end-1,:)).^2, 2));
% t_train = linspace(0, path_length/mean_speed, numel(x_train))';	% Uniform time spacing
t_train = [0; cumsum(path_length(:))/mean_speed];

fprintf(1, ['\n--- Path parameters ---\nPath length\t= %0.5gm\n',...
	'Time taken\t= %0.5gs\nMean speed\t= %0.5gm/s\n'], sum(path_length),...
	t_train(end), mean_speed);

W_train = W_actual(X_train', t_train)';
W_train = W_train + sigma_s*randn(size(X_train));

% -- Common, do not comment -- %
n_train = size(X_train, 1);

% -- Mean offset (steady wind overlay) -- %
W_mean = median(W_train, 1);
% W_train = W_train - ones(size(X_train, 1), 1)*W_mean;


%% TEST POINT DATA
x_test = linspace(min(X_train(:,1)), max(X_train(:,1)), 10)';
y_test = linspace(min(X_train(:,2)), max(X_train(:,2)), 10)';
z_test = linspace(min(X_train(:,3)), max(X_train(:,3)), 9)';

[x_grid, y_grid, z_grid] = meshgrid(x_test, y_test, z_test);
X_test = [x_grid(:), y_grid(:), z_grid(:)];


%% GP OPTIONS AND HYPERPARAMETER ESTIMATES

% Square-exponential spatial covariance
length_scale = 100;	% Spatial units
sigma_f = 1;
sigma_n = 0.1;		% NOTE: THIS COVERS TIME NOISE AS WELL!! Does it??

% Exponential temporal covariance
lag_scale = 10;		% Temporal units

% fh_covariance = covariance_sum(@square_expt, @d_square_expt, 
	
cov_funs = {@square_expt, @d_square_expt};
loghyper = log([length_scale, lag_scale, sigma_f, sigma_n]);

% cov_funs = {@cressie_ns5, @d_cressie_ns5};
% loghyper = [sigma_f, 200, 200, 1, sigma_n];

% Gaussian mean function
% NOTE: B values are standard deviation, which are squared to variance to
% prevent negative values.
mean_fun = @(X) ones([1, size(X, 1)]);
order_h = 1;

supercov_funs = {@square_expt};
superloghyper = [loghyper, ([W_mean(1), 1, W_mean(2), 1, W_mean(3), 1])];
supern_hyper = [4, 6];

% supercov_funs = {@neural_net_t};
% superloghyper = [ones(1,7), ([W_mean(1), 1, W_mean(2), 1, W_mean(3), 1])];
% supern_hyper = [7, 6];

% supercov_funs = {@cressie_ns1};
% superlohgyper = [[1, 50, 100, 0.1], ([W_mean(1), 1, W_mean(2), 1, W_mean(3), 1])];
% supern_hyper = [4, 6];

% supercov_funs = {@cressie_ns3};
% superloghyper = [2, 200, 200, 1, ([W_mean(1), 1, W_mean(2), 1, W_mean(3), 1])];
% supern_hyper = [4, 6];

% mean_fun = @(X) [ones([1, size(X, 1)]); X(:,1)'; X(:,2)'; X(:,3)'];
% superloghyper = [loghyper, repmat([0, ones(1,7)], 1, 3)];
% supern_hyper = [4, 24];

% Covariance sum version
cov_fun_handles = {@square_expt, @neural_net_cov};
ultran_hyper = [6, 5, 6];
ultraloghyper = [loghyper([1,1,1,2,3,4]), ones(1,5), ([W_mean(1), 1, W_mean(2), 1, W_mean(2), 1])];

% cov_fun_handles = {@neural_net_t, @neural_net_cov};
% ultran_hyper = [7, 5, 6];
% ultraloghyper = [ones(1,5), loghyper([2,4]), ones(1,5), ([W_mean(1), 1, W_mean(2), 1, W_mean(2), 1])];

ultracov_funs = {@(x1, t1, x2, t2, hyper) covariance_sum(x1, t1, x2, t2, hyper, cov_fun_handles, ultran_hyper)};


%% GP optimisation - GPt

opt = optimset('Display', 'iter', 'GradObj', 'on', 'TolX', 1e-3,...
	'MaxIter', 200);

% Spatio-temporal single covariance function
[fmin, nlml] = fminunc(@(hyper) ...
	GPt_likelihoodn(X_train, t_train, W_train(:,1), W_train(:,2), W_train(:,3), cov_funs, hyper), loghyper, opt);

disp('----- Results: Spatio-temporal squared exponential only -----');
fprintf(1, '\nMixed Solution: l = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n\n', ...
	exp(fmin(1)), exp(fmin(2)), exp(fmin(3)), exp(fmin(4)));


%% GP optimisation - GPt with gaussian mean function
opt = optimset('Display', 'iter', 'GradObj', 'off', 'TolX', 1e-5,...
	'MaxIter', 100);

[gfmin, gnlml] = fminunc(@(hyper) ...
	GPt_likelihoodn_gaussmean(X_train, t_train, W_train(:,1), W_train(:,2), ...
	W_train(:,3), supercov_funs, mean_fun, hyper, supern_hyper), superloghyper, opt);
% gfmin = positivify(gfmin, supern_hyper, order_h);

% Display results
disp('----- Results: Spatio-temporal squared exponential and gaussian mean function -----');
fprintf(1, '\nsq_expt function Solution: l = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(gfmin(1)), exp(gfmin(2)), exp(gfmin(3)), exp(gfmin(4)));
% fprintf(1, ['\nNeural_net_t function Solution: sigmaf = %0.5g, ',...
% 	'sigma0 = %0.5g, sigmax = %0.5g, sigmay = %0.5g, sigmaz = %0.5g,',...
% 	'sigma_t = %0.5g, sigma_n = %0.5g\n'], ...
%  	gfmin(1:5), exp(gfmin(6)), gfmin(7));
% fprintf(1, '\nCressie non-seperable 1 function Solution: sf = %0.5g, lt = %0.5g, l = %0.5g, sn = %0.5g\n', ...
% 	gfmin(1:4));
% fprintf(1, ['\nCressie non-seperable 5 function Solution: sf = %0.5g, ',...
% 	'lt = %0.5g, t = %0.5g, c = %0.5g, sn = %0.5g\n'], gfmin(1:5));


fprintf(1, 'Mean function:   '); fprintf(1, '%0.5g   ', gfmin(end-5:end));
fprintf(1, '\n\n');


%% GP optimisation - Multi-covariance GP + gaussian mean function
[ugfmin, ugnlml] = fminunc(@(hyper) ...
	GPt_likelihoodn_gaussmean(X_train, t_train, W_train(:,1), W_train(:,2), ...
	W_train(:,3), ultracov_funs, mean_fun, hyper, ultran_hyper), ultraloghyper, opt);
ugfmin = positivify(ugfmin, ultran_hyper, order_h);


% Display results
disp('----- Results: Spatio-temporal squared exponential, pure spatial square exponential and gaussian mean function -----');
% fprintf(1, '\nsq_expt function Solution: l = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
% 	exp(ugfmin(1:ultran_hyper(1))));
fprintf(1, '\nsq_exp function Solution: lx = %0.5g, ly = %0.5g, lz = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(ugfmin(1:ultran_hyper(1))));
fprintf(1, '\nNeural_net function Solution: sigmaf = %0.5g, sigma0 = %0.5g, sigmax = %0.5g, sigmay = %0.5g, sigmaz = %0.5g, \n', ...
 	(ugfmin((ultran_hyper(1)+1):(ultran_hyper(1)+ultran_hyper(2)))));
fprintf(1, 'Mean function:   ');
fprintf(1, '%0.5g   ', ugfmin(sum(ultran_hyper(1:end-1))+1:end));
fprintf(1, '\n\n');

%% <-> <-> <-> <-> <-> EDITING POSITION <-> <-> <-> <-> <-> %%

% [W_est, W_cov] = GPt_predict(X_train, t_train, W_train, ...
% 	X_test, t_train(1)*ones(size(X_test,1), 1), cov_funs{1}, fmin);

% [W_est, W_cov] = GPt_predict_gaussmean(X_train, t_train, W_train, ...
% 	X_test, t_train(1)*ones(size(X_test,1), 1), cov_funs{1}, mean_fun, gfmin);

% [W_est, W_cov] = GPt_predict_gaussmean(X_train, t_train, W_train, ...
% 	X_test, t_train(1)*ones(size(X_test,1), 1), ultracov_funs{1}, mean_fun, ugfmin, ultran_hyper);

% W_est = W_est + ones(size(W_est, 1), 1)*W_mean;

%% Plotting
% figure(1); clf;
% quiver3(x_grid, y_grid, z_grid, reshape(W_est(:,1), size(x_grid)), ...
% 	reshape(W_est(:,2), size(x_grid)), reshape(W_est(:,3), size(x_grid)))
% 
% hold on;
% W_grid = W_actual(X_test,t_train(1)*ones(size(X_test,1), 1) );
% quiver3(x_grid, y_grid, z_grid, reshape(W_grid(:,1), size(x_grid)), ...
% 	reshape(W_grid(:,2), size(x_grid)), reshape(W_grid(:,3), size(x_grid)))



%% Estimate videos
plot_actual = false;
plot_est = true;
remove_mean = false;
lt_index = 4;

figure(2); clf;
set(gcf, 'Position', [100, 100, 640, 480]);
daspect([1,1,1]); hold on;

W_grid = W_actual(X_test',t_train(1)*ones(size(X_test,1), 1) )';
max_wind = sqrt(max(sum(W_grid.^2, 2)));
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
[W_est, W_cov] = GPt_predictn_gaussmean(X_train, t_train, W_train, ...
	X_test, t_train(1)*ones(size(X_test,1), 1), ultracov_funs{1}, mean_fun, ugfmin, ultran_hyper);
max_cov = max(W_cov);


if remove_mean
W_est = (W_est - ones(size(W_est, 1), 1)*W_steady);	% REMOVE MEAN WIND AND RESCALE - DODGY
% W_est = (W_est - cos_profile3(X_test, cos_layer));	% REMOVE MEAN WIND
end

observe_color = square_exp(t_train(1), t_train, [ugfmin(lt_index), 0]);
h_train = scatter3(X_train(:,1), X_train(:,2), X_train(:,3), ...
	marker_size*observe_color, max_cov - max_cov*observe_color, 'filled');

est_scale = sqrt(max(sum(W_est.^2, 2)))/max_wind;
h_cone_est = coneplot(x_grid, y_grid, z_grid, ...
	reshape(W_est(:,1), size(x_grid)), ...
	reshape(W_est(:,2), size(x_grid)), ...
	reshape(W_est(:,3), size(x_grid)), ...
	x_grid, y_grid, z_grid, est_scale, reshape(W_cov, size(x_grid)));
set(h_cone_est, 'EdgeColor', 'none');

else
h_train = scatter3(X_train(:,1), X_train(:,2), X_train(:,3), ...
	marker_size, 'filled');
end


set(gca, 'ZDir', 'reverse', 'YDir', 'reverse'); axis equal; axis tight;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
M(1) = getframe(gcf);
M = repmat(M(1), 1, n_train);
% fprintf('\n\nFrame 1 of %g', n_train)

for i = 2:n_train
	
	if plot_est
	[W_est, W_cov] = GPt_predictn_gaussmean(X_train, t_train, W_train, ...
	X_test, t_train(i)*ones(size(X_test,1), 1), ultracov_funs{1}, mean_fun, ugfmin, ultran_hyper);

% 	max_est = sqrt(max(sum(W_grid.^2, 2)));
	
	if remove_mean
		W_est = (W_est - ones(size(W_est, 1), 1)*W_steady);	% REMOVE MEAN WIND - DODGY
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
	
	observe_color = square_exp(t_train(i), t_train, [ugfmin(lt_index), 0]);
	set(h_train, 'Cdata', max_cov - max_cov*observe_color, ...
		'SizeData', marker_size*observe_color);
	end
	
	if plot_actual
	W_grid = W_actual(X_test,t_train(i)*ones(size(X_test,1), 1) );

	if remove_mean
		W_grid = (W_grid - ones(size(W_grid, 1), 1)*W_steady);	% REMOVE MEAN WIND - DODGY
		% W_grid = (W_grid - cos_profile3(X_test, cos_layer));	% REMOVE MEAN WIND
	end
		
	delete(h_cone_actual);
	est_scale = sqrt(max(sum(W_grid.^2, 2)))/max_wind;
	h_cone_actual = coneplot(x_grid, y_grid, z_grid, ...
		reshape(W_grid(:,1), size(x_grid)), ...
		reshape(W_grid(:,2), size(x_grid)), ...
		reshape(W_grid(:,3), size(x_grid)), ...
		x_grid, y_grid, z_grid, est_scale);
	set(h_cone_actual, 'EdgeColor', 'none', 'facecolor', [.5 .5 .5]);
	end
	
	
	M(i) = getframe(gcf);
% 	fprintf('\rFrame %g of %g', i, n_train)
end

fprintf('\nAll frames complete; converting to avi...\n');
movie2avi(M, 'GPt_movies\OVERWRITE.avi', 'fps', 10, 'Compression','none')
fprintf('Finished\n');

%% SCRAP
% [W_est, W_cov] = GPt_predict_gaussmean(X_train, t_train, W_train, ...
% 	X_test, t_train(1)*ones(size(X_test,1), 1), cov_funs{1}, mean_fun, gfmin);


	% 	[W_est, W_cov] = GPt_predict(X_train, t_train, W_train, ...
	% 		X_test, t_train(i)*ones(size(X_test,1), 1), cov_funs{1}, fmin);
	
% 	[W_est, W_cov] = GPt_predict_gaussmean(X_train, t_train, W_train, ...
% 		X_test, t_train(i)*ones(size(X_test,1), 1), cov_funs{1}, mean_fun,
% 		gfmin);

% clear M W_cov W_est W_grid h_cone_actual h_cone_est h_train marker_size
% max_cov max_est max_wind observe_color