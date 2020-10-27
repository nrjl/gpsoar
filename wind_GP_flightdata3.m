%% Setup for running GP over flight wind data
clear variables
close all

%% GENERATE OR RETRIEVE FULL SET OF TRAINING DATA
% Get training data from log file. Log contains x and y coordinates stored
% as x and y vectors and corresponding u and v wind components

% ----- CUB DATA -----
% load Cub_005_data

% ----- MATLAB WIND DATA -----
% load wind; x=x(:); y=y(:); z=z(:); u=u(:); v=v(:);
% dstart = 1;
% dstop = numel(x);
% x = x(dstart:dstop);
% y = y(dstart:dstop);
% z = z(dstart:dstop);
% u = u(dstart:dstop);
% v = v(dstart:dstop);

% ----- THERMAL BUBBLE -----
k_thermal = 2;
x = -200:50:200;
[x, y, z] = meshgrid(x, x, k_thermal*x);
wind_function = @(x, y, z) torus_thermal(x, y, z, 3, 100, k_thermal, [0; 0; 0]);

[u, v, w] = wind_function(x, y, z);
x=x(:); y=y(:); z=z(:); u=u(:); v=v(:); w = w(:);

%% TRAINING DATA
% These are the training data points, i.e. the points that you have sampled
% and from which you wish to reconstruct the field.

% -- Uniform sampling from input data --

% Skip number of training set (ie use every nth data point)
% train_skip = 15;		
% 
% x_train = x(1:train_skip:end);
% y_train = y(1:train_skip:end);
% z_train = z(1:train_skip:end);
% 
% a_train = u(1:train_skip:end);
% b_train = v(1:train_skip:end);

% -- Generated path samples --

% Helical path - parametrised by t. Sinusoidal variation of radius.
ri = 0.25*max(abs(x(:)));
ro = 0.75*max(abs(x(:)));
h_base = max(z); h_top = min(z); delta_h = h_top-h_base;
n_loops = 6;
n_samples = 200;

h_radius = @(t, n, ri, ro) (ro-ri)/2*(cos(t/n)-1) + ro;
t = linspace(0, 2*pi*n_loops, n_samples);

x_train = h_radius(t, n_loops, ri, ro).*cos(t);
y_train = h_radius(t, n_loops, ri, ro).*sin(t);
z_train = h_base + delta_h/(2*pi*n_loops)*t;

[a_train, b_train, c_train] = wind_function(x_train, y_train, z_train);
a_train = sqrt(a_train.^2 + b_train.^2);
b_train = c_train;

% -- Common, do not comment -- %
n_train = numel(x_train);
train_set = [x_train(:), y_train(:), z_train(:)];
a_train = a_train(:); b_train = b_train(:);

%% TEST POINT DATA
% These are the points that you want to estimate the mean and covariance of
% the function for, i.e the points that you want to sample from the
% reconstructed function.

% Uniform test point field
nx_test = 10;
ny_test = 10;
nz_test = 10;

uv_coord = 1;

x_test = linspace(min(x), max(x), nx_test);
y_test = linspace(min(y), max(y), ny_test);
z_test = linspace(min(z), max(z), nz_test);
[x_test, y_test, z_test] = meshgrid(x_test, y_test, z_test);

test_set = [x_test(:), y_test(:), z_test(:)];

%% GP OPTIONS AND HYPERPARAMETER ESTIMATES
length_scale = 100;
sigma_f = 100;
sigma_n = 1;
cov_funs = {@square_exp, @d_square_exp};
loghyper = log([length_scale, sigma_f, sigma_n]);

%% Plot original field (which is effectively training data)
time_color = jet(n_train);

% Color scaled to wind values
colour_resolution = 100;
cmap = jet(colour_resolution+1);
getcolor = @(x, scale, res, map) map(round((scale(2) - double(x))/...
	(scale(2) - scale(1))*res+1) ,:);

fullcolor_a = getcolor(a_train, [min(a_train), max(a_train)], colour_resolution, cmap);
fullcolor_b = getcolor(b_train, [min(b_train), max(b_train)], colour_resolution, cmap);

figure(1); clf;
set(gcf, 'Name', 'Original data');
h_quiver1 = subplot(2,1,1);
quiver3(x, y, z, u, v, zeros(size(u)));
axis equal; axis tight; hold on;
title('Original wind vector field');
xlabel('X'); ylabel('Y'); zlabel('Z');

h_a1 = subplot(2,2,3); hold on;
% plot3(train_set(:,1), train_set(:,2), a_train, 'r+');
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	4*ones(n_train, 1), fullcolor_a, 'filled');
axis tight; view(3); title('a');
xlabel('X'); ylabel('Y'); zlabel('a');
	
h_b1 = subplot(2,2,4); hold on;
% plot3(train_set(:,1), train_set(:,2), b_train, 'r+');
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	4*ones(n_train, 1), fullcolor_b, 'filled');
axis tight; view(3); title('b');
xlabel('X'); ylabel('Y'); zlabel('b');

drawnow;
% return;


%% GP optimisation

opt = optimset('Display', 'iter', 'GradObj', 'on', 'TolX', 1e-7,...
	'MaxIter', 100);

% a optimise
[fmina, nlmla] = fminunc(...
	@(hyper) GP_likelihood(train_set, a_train, cov_funs, hyper), loghyper, opt);
fprintf(1, '\nSolution a: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fmina(1)), exp(fmina(2)), exp(fmina(3)));
[a_optimum, a_cov_optimum] = ...
	GP_predict(train_set, a_train, test_set, cov_funs{1}, fmina);

% b optimise
[fminb, nlmlb] = fminunc(...
	@(hyper) GP_likelihood(train_set, b_train, cov_funs, hyper), fmina, opt);

fprintf(1, '\nSolution b: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fminb(1)), exp(fminb(2)), exp(fminb(3)));
[b_optimum, b_cov_optimum] = ...
	GP_predict(train_set, b_train, test_set, cov_funs{1}, fminb);

% Group optimise
[fminab, nlmlab] = fminunc(@(hyper) ...
	GP_likelihoodn(train_set, a_train, b_train, cov_funs, hyper), loghyper, opt);
fprintf(1, '\nMixed Solution: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fminab(1)), exp(fminab(2)), exp(fminab(3)));
[a_mix, a_cov_mix] = ...
	GP_predict(train_set, a_train, test_set, cov_funs{1}, fminab);
[b_mix, b_cov_mix] = ...
	GP_predict(train_set, b_train, test_set, cov_funs{1}, fminab);

%% Plot setup

a_covariance = reshape(a_cov_optimum, size(x_test));
b_covariance = reshape(b_cov_optimum, size(x_test));

a_star = reshape(a_optimum, size(x_test));
b_star = reshape(b_optimum, size(x_test));

a_std = sqrt(a_covariance);
b_std = sqrt(b_covariance);

correct = @(a, amin, amax, bmin, bmax) ...
	(a-amin)/(amax-amin)*(bmax-bmin) + bmin;

a_colourcorrect = correct(a_std, min(a_std(:)), max(a_std(:)), ...
	min(a_star(:)), max(a_star(:)));
b_colourcorrect = correct(b_std, min(b_std(:)), max(b_std(:)), ...
	min(b_star(:)), max(b_star(:)));


if ~uv_coord
	u_est = a_star.*cos(b_star);
	v_est = -a_star.*sin(b_star);
else
	u_est = a_star;
	v_est = b_star;
end

arrowscale = max(sqrt(u_est(:).^2 + v_est(:).^2))/...
	max(sqrt(u(:).^2 + v(:).^2));

% Get overall boundaries from training and test data
lims_a = [min([min(a_optimum), min(a_train)]), max([max(a_optimum), max(a_train)])];
lims_b = [min([min(b_optimum), min(b_train)]), max([max(b_optimum), max(b_train)])]; 

% Get colors of test points
fullcolor_atest = getcolor(a_optimum, lims_a, colour_resolution, cmap);
fullcolor_btest = getcolor(b_optimum, lims_b, colour_resolution, cmap);

% Get colors of training points with matched colour limits
fullcolor_amatch = getcolor(a_train, lims_a, colour_resolution, cmap);
fullcolor_bmatch = getcolor(b_train, lims_b, colour_resolution, cmap);

%%
figure(2); clf;
set(gcf, 'Name', 'Separate solutions');
h_quiver2 = subplot(2,1,1);
% contour(x_test, y_test, a_std);

% quiver3(x_test, y_test, z_test, u_est, v_est, zeros(size(u_est)),
% arrowscale);
daspect([1,1,1]);
h_cone2 = coneplot(x_test, y_test, z_test, u_est, v_est, ...
	zeros(size(u_est)), x_test, y_test, z_test, a_std);
set(h_cone2, 'EdgeColor', 'none')
axis tight; hold on;
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');

h_a2 = subplot(2,2,3);
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_atest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_amatch, 'filled');
axis tight;
title('GP Estimated - Variable a');
xlabel('X'); ylabel('Y'); zlabel('Z');

h_b2 = subplot(2,2,4); 
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_btest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_bmatch, 'filled');
axis tight;
title('GP Estimated - Variable b');
xlabel('X'); ylabel('Y'); zlabel('Z');

%% Plot mixed solution

mix_covariance = reshape(a_cov_mix, size(x_test));

a_mix = reshape(a_mix, size(x_test));
b_mix = reshape(b_mix, size(x_test));

mix_std = sqrt(mix_covariance);

correct = @(a, amin, amax, bmin, bmax) ...
	(a-amin)/(amax-amin)*(bmax-bmin) + bmin;

mix_colourcorrect = correct(mix_std, min(mix_std(:)), max(mix_std(:)), ...
	min(a_star(:)), max(a_star(:)));

yellow = [1, 1, 0.85];
red = [1, 0.9, 0.9];

if ~uv_coord
	u_mix = a_mix.*cos(b_mix);
	v_mix = -a_mix.*sin(b_mix);
else
	u_mix = a_mix;
	v_mix = b_mix;
end

arrowscale = max(sqrt(u_mix(:).^2 + v_mix(:).^2))/...
	max(sqrt(u(:).^2 + v(:).^2));

% Get overall boundaries from training and test data
lims_amix = [min([min(a_mix(:)), min(a_train)]), max([max(a_mix(:)), max(a_train)])];
lims_bmix = [min([min(b_mix(:)), min(b_train)]), max([max(b_mix(:)), max(b_train)])]; 

% Get colors of test points
fullcolor_amixtest = getcolor(a_mix, lims_amix, colour_resolution, cmap);
fullcolor_bmixtest = getcolor(b_mix, lims_bmix, colour_resolution, cmap);

% Get colors of training points with matched colour limits
fullcolor_amixmatch = getcolor(a_train, lims_amix, colour_resolution, cmap);
fullcolor_bmixmatch = getcolor(b_train, lims_bmix, colour_resolution, cmap);

%%
figure(3); clf;
set(gcf, 'Name', 'Mixed solution');

h_quiver3 = subplot(2,1,1);
daspect([1,1,1]);
h_cone3 = coneplot(x_test, y_test, z_test, u_mix, v_mix, ...
	zeros(size(u_mix)), x_test, y_test, z_test, mix_std);
set(h_cone3, 'EdgeColor', 'none')
axis tight; hold on;
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');

h_a3 = subplot(2,2,3);
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_amixtest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_amixmatch, 'filled');
axis tight;
title('GP Estimated - Variable a');
xlabel('X'); ylabel('Y'); zlabel('Z');

h_b3 = subplot(2,2,4); 
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_bmixtest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_bmixmatch, 'filled');
axis tight;
title('GP Estimated - Variable b');
xlabel('X'); ylabel('Y'); zlabel('Z');



%%

linkaxes([h_quiver1, h_quiver2, h_quiver3]);

hlinka = linkprop([h_a3, h_a2, h_a1], {'Xlim', 'Ylim', 'Zlim', 'View'});
hlinkb = linkprop([h_b3, h_b2, h_b1], {'Xlim', 'Ylim', 'Zlim', 'View'});


%% List results table

% Calculate estimation at all training points
[a_full_est] = ...
	GP_predict(train_set, a_train, [x(:), y(:)], cov_funs{1}, fmina);
[b_full_est] = ...
	GP_predict(train_set, b_train, [x(:), y(:)], cov_funs{1}, fminb);
[a_full_mix] = ...
	GP_predict(train_set, a_train, [x(:), y(:)], cov_funs{1}, fminab);
[b_full_mix] = ...
	GP_predict(train_set, b_train, [x(:), y(:)], cov_funs{1}, fminab);

err_a = sum((a_full_est(:) - u(:)).^2);
err_b = sum((b_full_est(:) - v(:)).^2);
err_ab = sum((a_full_mix(:) - u(:)).^2) + sum((b_full_mix(:) - v(:)).^2);

%%
fprintf('\n----- GP Training results -----\n');
fprintf(['Method\t\t| length\t| sigma_f\t| sigma_n\t| nlml\t\t|| ', ...
    'Sum Squared Error\n']);
disp(repmat('-', 1, 80));
fprintf('a (alone)\t| %6.3f\t| %6.5f\t| %6.5f\t| %6.5f\t|| %6.5f\n', ...
	exp(fmina), nlmla, err_a);
fprintf('b (alone)\t| %6.3f\t| %6.5f\t| %6.5f\t| %6.5f\t|| %6.5f\n', ...
	exp(fminb), nlmlb, err_b);
fprintf('Grouped\t\t| %6.3f\t| %6.5f\t| %6.5f\t| %6.5f\t|| %6.5f\n', ...
	exp(fminab), nlmlab, err_ab);


%% SCRAP
% h_quiver3 = subplot(2,1,1);
% contour(x_test, y_test, mix_std);
% axis equal; axis tight; hold on;
% quiver(x_test, y_test, u_est, v_est, arrowscale);
% title('Estimated wind vector field');
% xlabel('X'); ylabel('Y');
% plot(h_quiver3, train_set(:,1), train_set(:,2), 'r+');
%
% h_a3 = subplot(2,2,3);
% surf(x_test, y_test, a_mix); hold on;
% surf(x_test, y_test, a_mix + 2*mix_std, a_colourcorrect, ...
% 	'FaceAlpha', 0.5, 'EdgeColor', 'none');
% surf(x_test, y_test, a_mix - 2*mix_std, a_colourcorrect, ...
% 	'FaceAlpha', 0.5, 'EdgeColor', 'none');
% title('GP Estimated - Variable a');
% xlabel('X'); ylabel('Y'); zlabel('a');
% plot3(h_a3, train_set(:,1), train_set(:,2), a_train, 'r+');
% axis tight;
% 
% h_b3 = subplot(2,2,4); 
% surf(x_test, y_test, b_mix); hold on;
% surf(x_test, y_test, b_mix + 2*mix_std, b_colourcorrect, ...
% 	'FaceAlpha', 0.5, 'EdgeColor', 'none');
% surf(x_test, y_test, b_mix - 2*mix_std, b_colourcorrect, ...
% 	'FaceAlpha', 0.5, 'EdgeColor', 'none');
% title('GP Estimated - Variable b');
% xlabel('X'); ylabel('Y'); zlabel('b');
% plot3(h_b3, train_set(:,1), train_set(:,2), b_train, 'r+');
% axis tight;
