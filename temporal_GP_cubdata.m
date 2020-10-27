%% Setup for running GP over flight wind data
clear variables
close all

plot_on = 1;
%% GENERATE OR RETRIEVE FULL SET OF TRAINING DATA

load cub_005_datat.mat

%% TRAINING DATA
% These are the training data points, i.e. the points that you have sampled
% and from which you wish to reconstruct the field.

pskip = 20;
istart = 1000;
istop = 6000;

X_train = [x(istart:pskip:istop), y(istart:pskip:istop), z(istart:pskip:istop)];
t_train = t(istart:pskip:istop)-t(1);

path_length = sqrt(sum((X_train(2:end,:) - X_train(1:end-1,:)).^2, 2));

mean_speed = sum(path_length)/(t_train(end)-t_train(1));

fprintf(1, ['\n--- Path parameters ---\nPath length\t= %0.5gm\n',...
	'Time taken\t= %0.5gs\nMean speed\t= %0.5gm/s\n'], sum(path_length),...
	t_train(end), mean_speed);

W_train = [u(istart:pskip:istop), v(istart:pskip:istop)];

% -- Common, do not comment -- %
n_train = size(X_train, 1);

% -- Mean offset (steady wind overlay) -- %
W_mean = median(W_train, 1);


%% TEST POINT DATA
x_test = linspace(min(X_train(:,1)), max(X_train(:,1)), 10)';
y_test = linspace(min(X_train(:,2)), max(X_train(:,2)), 10)';
z_test = linspace(min(X_train(:,3)), max(X_train(:,3)), 9)';

[x_grid, y_grid, z_grid] = meshgrid(x_test, y_test, z_test);
X_test = [x_grid(:), y_grid(:), z_grid(:)];


%% GP OPTIONS AND HYPERPARAMETER ESTIMATES

% Square-exponential spatial covariance
length_scale = 80;	% Spatial units
sigma_f = 1;
sigma_n = 0.1;		% NOTE: THIS COVERS TIME NOISE AS WELL!! Does it??

% Exponential temporal covariance
lag_scale = 30;		% Temporal units

cov_funs = {@square_expt, @d_square_expt};
loghyper = log([length_scale, length_scale, length_scale, lag_scale, sigma_f, sigma_n]);

% cov_funs = {@cressie_ns5, @d_cressie_ns5};
% loghyper = [sigma_f, 200, 200, 1, sigma_n];

% Gaussian mean function
% NOTE: B values are standard deviation, which are squared to variance to
% prevent negative values.
mean_fun = @(X) ones([1, size(X, 1)]);
order_h = 1;

supercov_funs = {@square_expt};
superloghyper = [loghyper, W_mean(1), 1, W_mean(2), 1];
supern_hyper = [6, 4];

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
ultran_hyper = [6, 5, 4];
ultraloghyper = [loghyper, ones(1,5), W_mean(1), 1, W_mean(2), 1];

% cov_fun_handles = {@neural_net_t, @neural_net_cov};
% ultran_hyper = [7, 5, 6];
% ultraloghyper = [ones(1,5), loghyper([2,4]), ones(1,5), ([W_mean(1), 1, W_mean(2), 1, W_mean(2), 1])];

ultracov_funs = {@(x1, t1, x2, t2, hyper) covariance_sum(x1, t1, x2, t2, hyper, cov_fun_handles, ultran_hyper)};


%% GP optimisation - GPt

opt = optimset('Display', 'iter', 'GradObj', 'on', 'TolX', 1e-3,...
	'MaxIter', 200);

% Spatio-temporal single covariance function
[fmin, nlml] = fminunc(@(hyper) ...
	GPt_likelihoodn(X_train, t_train, W_train(:,1), W_train(:,2), cov_funs, hyper), loghyper, opt);

disp('----- Results: Spatio-temporal squared exponential only -----');
fprintf(1, '\nMixed Solution: l = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n\n', ...
	exp(fmin));


%% GP optimisation - GPt with gaussian mean function
opt = optimset('Display', 'iter', 'GradObj', 'off', 'TolX', 1e-5,...
	'MaxIter', 100);

[gfmin, gnlml] = fminunc(@(hyper) ...
	GPt_likelihoodn_gaussmean(X_train, t_train, W_train(:,1), W_train(:,2), ...
	supercov_funs, mean_fun, hyper, supern_hyper), superloghyper, opt);
% gfmin = positivify(gfmin, supern_hyper, order_h);

% Display results
disp('----- Results: Spatio-temporal squared exponential and gaussian mean function -----');
% fprintf(1, '\nsq_expt function Solution: l = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
% 	exp(gfmin(1)), exp(gfmin(2)), exp(gfmin(3)), exp(gfmin(4)));
fprintf(1, '\nsq_exp function Solution: lx = %0.5g, ly = %0.5g, lz = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(ugfmin(1:supern_hyper(1))));
% fprintf(1, ['\nNeural_net_t function Solution: sigmaf = %0.5g, ',...
% 	'sigma0 = %0.5g, sigmax = %0.5g, sigmay = %0.5g, sigmaz = %0.5g,',...
% 	'sigma_t = %0.5g, sigma_n = %0.5g\n'], ...
%  	gfmin(1:5), exp(gfmin(6)), gfmin(7));
% fprintf(1, '\nCressie non-seperable 1 function Solution: sf = %0.5g, lt = %0.5g, l = %0.5g, sn = %0.5g\n', ...
% 	gfmin(1:4));
% fprintf(1, ['\nCressie non-seperable 5 function Solution: sf = %0.5g, ',...
% 	'lt = %0.5g, t = %0.5g, c = %0.5g, sn = %0.5g\n'], gfmin(1:5));


fprintf(1, 'Mean function:   ');
fprintf(1, '%0.5g   ', gfmin(sum(supern_hyper(1:end-1))+1:end));
fprintf(1, '\n\n');

%% GP optimisation - Multi-covariance GP + gaussian mean function
[ugfmin, ugnlml] = fminunc(@(hyper) ...
	GPt_likelihoodn_gaussmean(X_train, t_train, W_train(:,1), W_train(:,2), ...
	ultracov_funs, mean_fun, hyper, ultran_hyper), ultraloghyper, opt);


% Display results
disp('----- Results: Spatio-temporal squared exponential, pure spatial square exponential and gaussian mean function -----');
% fprintf(1, '\nsq_expt function Solution: l = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
% 	exp(ugfmin(1:ultran_hyper(1))));
fprintf(1, '\nsq_exp function Solution: lx = %0.5g, ly = %0.5g, lz = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(ugfmin(1:ultran_hyper(1))));
fprintf(1, '\nNeural_net function Solution: sigmaf = %0.5g, sigma0 = %0.5g, sigmax = %0.5g, sigmay = %0.5g, sigmaz = %0.5g\n', ...
 	(ugfmin((ultran_hyper(1)+1):(ultran_hyper(1)+ultran_hyper(2)))));
fprintf(1, 'Mean function:   ');
fprintf(1, '%0.5g   ', ugfmin(sum(ultran_hyper(1:end-1))+1:end));
fprintf(1, '\n\n');