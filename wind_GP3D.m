%% Setup for running GP over flight wind data
clear variables
close all
addpath gp_tools wind_functions

plot_on = 1;
%% GENERATE OR RETRIEVE FULL SET OF TRAINING DATA
% Get training data from log file. Log contains x and y coordinates stored
% as x and y vectors and corresponding u and v wind components

% ----- CUB DATA -----
% load Cub_005_data
% w = zeros(size(u));

% ----- MATLAB WIND DATA -----
% load wind; data_skip = 3;
% x_full = x; y_full = y; z_full = z;
% u_full = u; v_full = v; w_full = w;
% 
% x_grid = x(1:data_skip:end, 1:data_skip:end, 1:data_skip:end);
% y_grid = y(1:data_skip:end, 1:data_skip:end, 1:data_skip:end);
% z_grid = z(1:data_skip:end, 1:data_skip:end, 1:data_skip:end);
% 
% u_grid = u(1:data_skip:end, 1:data_skip:end, 1:data_skip:end);
% v_grid = v(1:data_skip:end, 1:data_skip:end, 1:data_skip:end);
% w_grid = w(1:data_skip:end, 1:data_skip:end, 1:data_skip:end);
% 
% x = x_grid(:); y = y_grid(:); z = z_grid(:);
% u = u_grid(:); v = v_grid(:); w = w_grid(:);
% dstart = 1;
% dstop = numel(x);
% x = x(dstart:dstop); y = y(dstart:dstop); z = z(dstart:dstop);
% u = u(dstart:dstop); v = v(dstart:dstop); w = w(dstart:dstop);

% ----- SINGLE THERMAL BUBBLE -----
% R_thermal = 100;
% k_thermal = 2;
% x = linspace(-2*R_thermal, 2*R_thermal, 9);
% [x_grid, y_grid, z_grid] = meshgrid(x, x, k_thermal/2*x);
% wind_function = @(x, y, z) torus_thermal(x, y, z, -3, R_thermal, ...
% k_thermal, [50; 0; 0]);

% ----- THERMAL FIELD ----- 
t1 = @(x, y, z) torus_thermal(x, y, z, -5, 50, 2, [50;0;50]);
t2 = @(x, y, z) torus_thermal(x, y, z, -3, 100, 2, [0;300;100]);
gr =  @(x, y, z) pohlhausen(z, 0, 2, 200, 0, -45);
% gr = @(x, y, z) sine_wind(x, y, z, 2, 500, 500);

x = linspace(-200, 200, 7);
y = [linspace(-100, 100, 6), 150, linspace(200, 400, 6)];
z = linspace(0, 250, 5);
[x_grid, y_grid, z_grid] = meshgrid(x, y, z);
wind_function = @(x, y, z) multi_field(x, y, z, t1, t2, gr);

% ----- SINUSOIDAL FIELD -----
% % V_sine = 3;
% % V_sine = @(x, y, z) pohlhausen(z, 110, 12, 50, 0);
% V_sine = @(x, y, z) linear_profile(z, 120, 170, 0, 15);
% L_xy = 500;
% L_xz = 500;
% x = linspace(0, 500, 10);
% y = linspace(-50, 50, 10);
% z = linspace(100, 200, 10);
% 
% [x_grid, y_grid, z_grid] = meshgrid(x, y, z);
% wind_function = @(x, y, z) sine_wind(x, y, z, V_sine, L_xy, L_xz);

% ----- Common for wind fields ----- %
[u_grid, v_grid, w_grid] = wind_function(x_grid, y_grid, z_grid);
x = x_grid(:); y = y_grid(:); z = z_grid(:);
u = u_grid(:); v = v_grid(:); w = w_grid(:);


%% TRAINING DATA
% These are the training data points, i.e. the points that you have sampled
% and from which you wish to reconstruct the field.

% -- Uniform sampling from input data --

% Skip number of training set (ie use every nth data point)
% train_skip = 10;		
% 
% x_train = x(1:train_skip:end);
% y_train = y(1:train_skip:end);
% z_train = z(1:train_skip:end);
% 
% a_train = u(1:train_skip:end);
% b_train = v(1:train_skip:end);
% c_train = w(1:train_skip:end);

% --- Pre generated path with interpolation to nearest point --- %
% WRITE IT YOURSELF %

% --- Generated path samples --- %

% -- Curvy path for Matlab wind data -- %
% n_cycles = 1.5;
% n_samples = 200;
% 
% % Data coverage of path sweep
% x_cover = 1.0;		y_cover = 0.8;		z_cover = 0.7;
% minx = min(x);		delx = max(x)-minx;
% miny = min(y);		dely = max(y)-miny;
% minz = min(z);		delz = max(z)-minz;
% 
% t = linspace(0, 2*n_cycles*pi, n_samples);
% x_train = minx + 0.5*(1-x_cover)*delx + x_cover*delx*t/(2*pi*n_cycles);
% y_train = miny + dely/2 - y_cover*dely/2*sin(t);
% z_train = minz + delz/2 + z_cover*delz/2*cos(t);


% -- Helical path -- %
% ri = 0.30*max(x);
% ro = 0.50*max(x);
% zb = max(z); zt = min(z);
% n_cycles = 4;
% n_samples = 200;

% [x_train, y_train, z_train] = ...
% 	helical_path(ri, ro, zb, zt, n_cycles, n_samples);

% -- Path through thermal field -- %
% Hand coded so won't automatically generate for a random thermal field
t0 = linspace(0, -pi, 10);
x0 = [linspace(-150, 160, 20), 110+50*cos(t0)]; 
y0 = [linspace(400, 0, 20), 50*sin(t0)];
z0 = linspace(20, 0, numel(x0));
[x1, y1, z1] = helical_path(20, 60, 0, 80, 3, 75);
x2 = linspace(60, 120, 10);
y2 = linspace(0, 300, 10);
z2 = linspace(80, 40, 10);
[x3, y3, z3] = helical_path(50, 120, 40, 250, 4, 100);
x_train = [x0, x1, x2, x3];
y_train = [y0, y1, y2, y3+300];
z_train = [z0, z1, z2, z3];

% Lookup or calculate wind values
try functions(wind_function);	% If wind function exists, then evaluate
	[a_train, b_train, c_train] = wind_function(x_train, y_train, z_train);

catch CE						% Otherwise lookup with nearest neighbour
	x_train = interp1(squeeze(x_full(1,:,1)), squeeze(x_full(1,:,1)), x_train, 'nearest');
	y_train = interp1(squeeze(y_full(:,1,1)), squeeze(y_full(:,1,1)), y_train, 'nearest');
	z_train = interp1(squeeze(z_full(1,1,:)), squeeze(z_full(1,1,:)), z_train, 'nearest');
	
	a_train = interp3(x_full, y_full, z_full, u_full, x_train, y_train, z_train, 'nearest');
	b_train = interp3(x_full, y_full, z_full, v_full, x_train, y_train, z_train, 'nearest');
	c_train = interp3(x_full, y_full, z_full, w_full, x_train, y_train, z_train, 'nearest');
end

% If full pose information (phi, theta, psi) is available then should
% convert to alpha, beta, V measurements, then add noise to measurements
% and restore, but currently just add noise directly
a_train = a_train + 0.1*randn(size(a_train));
b_train = b_train + 0.1*randn(size(b_train));
c_train = c_train + 0.1*randn(size(c_train));

% -- Common, do not comment -- %
n_train = numel(x_train);
train_set = [x_train(:), y_train(:), z_train(:)];
a_train = a_train(:); b_train = b_train(:);  c_train = c_train(:);

% -- Mean offset (steady wind overlay) -- %
mean_a = median(a_train); a_train = a_train - mean_a;
mean_b = median(b_train); b_train = b_train - mean_b;
mean_c = median(c_train); c_train = c_train - mean_c;

%% TEST POINT DATA
% These are the points that you want to estimate the mean and covariance of
% the function for, i.e the points that you want to sample from the
% reconstructed function.

% Uniform test point field
% nx_test = 9;
% ny_test = 9;
% nz_test = 9;
% 
% 
% x_test = linspace(min(x), max(x), nx_test);
% y_test = linspace(min(y), max(y), ny_test);
% z_test = linspace(min(z), max(z), nz_test);
% [x_test, y_test, z_test] = meshgrid(x_test, y_test, z_test);

% Same test set as grid set
x_test = x_grid; y_test = y_grid; z_test = z_grid;

test_set = [x_test(:), y_test(:), z_test(:)];
uv_coord = 1;
%% Plot original field (which is effectively training data)
time_color = jet(n_train);

% Color scaled to wind values
colour_resolution = 100;
cmap = jet(colour_resolution+1);
getcolor = @(x, scale, res, map) map(round((scale(2) - double(x))/...
	(scale(2) - scale(1))*res+1) ,:);

color_lim_u = [min([a_train; u]), max([a_train;u])];
color_lim_v = [min([b_train; v]), max([b_train;v])];
color_lim_w = [min([c_train; w]), max([c_train;w])];

fullcolor_a = getcolor(a_train, color_lim_u, colour_resolution, cmap);
fullcolor_b = getcolor(b_train, color_lim_v, colour_resolution, cmap);
fullcolor_c = getcolor(c_train, color_lim_w, colour_resolution, cmap);

fullcolor_u = getcolor(u, color_lim_u, colour_resolution, cmap);
fullcolor_v = getcolor(v, color_lim_v, colour_resolution, cmap);
fullcolor_w = getcolor(w, color_lim_w, colour_resolution, cmap);

if plot_on
figure(1); clf;
set(gcf, 'Name', 'Original data');

h_quiver1 = subplot(2,2,1);
% quiver3(x, y, z, u, v, w);
coneplot(x_grid, y_grid, z_grid, u_grid, v_grid, w_grid, ...
	x_grid, y_grid, z_grid);
axis equal; axis tight; hold on;
plot3(x_train, y_train, z_train, 'r+');
view(3);
title('Original wind vector field');
xlabel('X'); ylabel('Y'); zlabel('Z');

h_a1 = subplot(2,2,2); hold on;
scatter3(x, y, z, 4, fullcolor_u, 'filled');
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	4*ones(n_train, 1), fullcolor_a, 'filled');
axis tight; view(3); title('a');
xlabel('X'); ylabel('Y'); zlabel('a');
	
h_b1 = subplot(2,2,3); hold on;
scatter3(x, y, z, 4, fullcolor_v, 'filled');
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	4*ones(n_train, 1), fullcolor_b, 'filled');
axis tight; view(3); title('b');
xlabel('X'); ylabel('Y'); zlabel('b');

h_c1 = subplot(2,2,4); hold on;
scatter3(x, y, z, 4, fullcolor_w, 'filled');
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	4*ones(n_train, 1), fullcolor_c, 'filled');
axis tight; view(3); title('c');
xlabel('X'); ylabel('Y'); zlabel('c');

drawnow;
end
% return;

%% GP OPTIONS AND HYPERPARAMETER ESTIMATES
length_scale = 50; %(max(x(:)) - min(x(:)))/4;
sigma_f = 1;
sigma_n = 0.1;
cov_funs = {@square_exp, @d_square_exp};
loghyper = log([length_scale, sigma_f, sigma_n]);

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
a_optimum = a_optimum + mean_a;

% b optimise
[fminb, nlmlb] = fminunc(...
	@(hyper) GP_likelihood(train_set, b_train, cov_funs, hyper), loghyper, opt);

fprintf(1, '\nSolution b: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fminb(1)), exp(fminb(2)), exp(fminb(3)));
[b_optimum, b_cov_optimum] = ...
	GP_predict(train_set, b_train, test_set, cov_funs{1}, fminb);
b_optimum = b_optimum + mean_b;

% c optimise
[fminc, nlmlc] = fminunc(...
	@(hyper) GP_likelihood(train_set, c_train, cov_funs, hyper), loghyper, opt);

fprintf(1, '\nSolution c: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fminc(1)), exp(fminc(2)), exp(fminc(3)));
[c_optimum, c_cov_optimum] = ...
	GP_predict(train_set, c_train, test_set, cov_funs{1}, fminc);
c_optimum = c_optimum + mean_c;

% Group optimise
[fminab, nlmlab] = fminunc(@(hyper) ...
	GP_likelihoodn(train_set, a_train, b_train, c_train, cov_funs, hyper), loghyper, opt);
fprintf(1, '\nMixed Solution: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fminab(1)), exp(fminab(2)), exp(fminab(3)));
[abc_mix, cov_mix] = GP_predict(train_set, [a_train,b_train,c_train], ...
	test_set, cov_funs{1}, fminab);

a_mix = abc_mix(:,1) + mean_a;
b_mix = abc_mix(:,2) + mean_b;
c_mix = abc_mix(:,3) + mean_c;


%% Plot setup

a_covariance = reshape(a_cov_optimum, size(x_test));
b_covariance = reshape(b_cov_optimum, size(x_test));
c_covariance = reshape(c_cov_optimum, size(x_test));

a_star = reshape(a_optimum, size(x_test));
b_star = reshape(b_optimum, size(x_test));
c_star = reshape(c_optimum, size(x_test));

a_std = sqrt(a_covariance);
b_std = sqrt(b_covariance);
c_std = sqrt(c_covariance);

correct = @(a, amin, amax, bmin, bmax) ...
	(a-amin)/(amax-amin)*(bmax-bmin) + bmin;

a_colourcorrect = correct(a_std, min(a_std(:)), max(a_std(:)), ...
	min(a_star(:)), max(a_star(:)));
b_colourcorrect = correct(b_std, min(b_std(:)), max(b_std(:)), ...
	min(b_star(:)), max(b_star(:)));
c_colourcorrect = correct(c_std, min(c_std(:)), max(c_std(:)), ...
	min(c_star(:)), max(c_star(:)));


if ~uv_coord
	u_est = a_star.*cos(b_star);
	v_est = -a_star.*sin(b_star);
	w_est = c_star;
else
	u_est = a_star;
	v_est = b_star;
	w_est = c_star;
end

arrowscale = max(sqrt(u_est(:).^2 + v_est(:).^2 + w_est(:).^2))/...
	max(sqrt(u(:).^2 + v(:).^2 + w(:).^2));

% Get overall boundaries from training and test data
lims_a = [min([min(a_optimum), min(a_train)]), max([max(a_optimum), max(a_train)])];
lims_b = [min([min(b_optimum), min(b_train)]), max([max(b_optimum), max(b_train)])]; 
lims_c = [min([min(c_optimum), min(c_train)]), max([max(c_optimum), max(c_train)])]; 

% Get colors of test points
fullcolor_atest = getcolor(a_optimum, lims_a, colour_resolution, cmap);
fullcolor_btest = getcolor(b_optimum, lims_b, colour_resolution, cmap);
fullcolor_ctest = getcolor(c_optimum, lims_c, colour_resolution, cmap);

% Get colors of training points with matched colour limits
fullcolor_amatch = getcolor(a_train, lims_a, colour_resolution, cmap);
fullcolor_bmatch = getcolor(b_train, lims_b, colour_resolution, cmap);
fullcolor_cmatch = getcolor(c_train, lims_c, colour_resolution, cmap);

%% Plot separate results
zerob = zeros(size(u_est));

if plot_on
figure(2); clf;
set(gcf, 'Name', 'Seperate solutions');
h_quiver2 = subplot(2,2,1);
quiver3(x_test, y_test, z_test, u_est, v_est, w_est);


axis equal; axis tight; hold on;
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');
plot3(h_quiver2, train_set(:,1), train_set(:,2), train_set(:,3), 'r+');

h_a2 = subplot(2,2,2);
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_atest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_amatch, 'filled');
axis tight;
title('GP Estimated - Variable a');
xlabel('X'); ylabel('Y'); zlabel('Z');

h_b2 = subplot(2,2,3); 
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_btest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_bmatch, 'filled');
axis tight;
title('GP Estimated - Variable b');
xlabel('X'); ylabel('Y'); zlabel('Z');

h_c2 = subplot(2,2,4); 
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_ctest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_cmatch, 'filled');
axis tight;
title('GP Estimated - Variable c');
xlabel('X'); ylabel('Y'); zlabel('Z');
end

%% Plot mixed solution

cov_mix = reshape(cov_mix, size(x_test));

a_mix = reshape(a_mix, size(x_test));
b_mix = reshape(b_mix, size(x_test));
c_mix = reshape(c_mix, size(x_test));

mix_std = sqrt(cov_mix);

correct = @(a, amin, amax, bmin, bmax) ...
	(a-amin)/(amax-amin)*(bmax-bmin) + bmin;

mix_colourcorrect = correct(mix_std, min(mix_std(:)), max(mix_std(:)), ...
	min(a_mix(:)), max(a_mix(:)));

yellow = [1, 1, 0.85];
red = [1, 0.9, 0.9];

if ~uv_coord
	u_mix = a_mix.*cos(b_mix);
	v_mix = -a_mix.*sin(b_mix);
	w_mix = c_mix;
else
	u_mix = a_mix;
	v_mix = b_mix;
	w_mix = c_mix;
end

arrowscale_mix = max(sqrt(u_mix(:).^2 + v_mix(:).^2 + w_mix(:).^2))/...
	max(sqrt(u(:).^2 + v(:).^2 + w(:).^2));

% Get overall boundaries from training and test data
lims_amix = [min([min(a_mix(:)), min(a_train)]), max([max(a_mix(:)), max(a_train)])];
lims_bmix = [min([min(b_mix(:)), min(b_train)]), max([max(b_mix(:)), max(b_train)])]; 
lims_cmix = [min([min(c_mix(:)), min(c_train)]), max([max(c_mix(:)), max(c_train)])]; 

% Get colors of test points
fullcolor_amixtest = getcolor(a_mix, lims_amix, colour_resolution, cmap);
fullcolor_bmixtest = getcolor(b_mix, lims_bmix, colour_resolution, cmap);
fullcolor_cmixtest = getcolor(c_mix, lims_cmix, colour_resolution, cmap);

% Get colors of training points with matched colour limits
fullcolor_amixmatch = getcolor(a_train, lims_amix, colour_resolution, cmap);
fullcolor_bmixmatch = getcolor(b_train, lims_bmix, colour_resolution, cmap);
fullcolor_cmixmatch = getcolor(c_train, lims_cmix, colour_resolution, cmap);

%%
if plot_on
figure(3); clf;
set(gcf, 'Name', 'Common hyperparameters solution');
h_quiver3 = subplot(2,2,1);
quiver3(x_test, y_test, z_test, u_mix, v_mix, w_mix);
axis equal; axis tight; hold on;
title('Estimated wind vector field');
xlabel('X'); ylabel('Y'); zlabel('Z');
plot3(h_quiver3, train_set(:,1), train_set(:,2), train_set(:,3), 'r+');

h_a3 = subplot(2,2,2);
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_bmixtest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_amixmatch, 'filled');
axis tight;
title('GP Estimated - Variable a');
xlabel('X'); ylabel('Y'); zlabel('Z');

h_b3 = subplot(2,2,3); 
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_bmixtest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_bmixmatch, 'filled');
axis tight;
title('GP Estimated - Variable b');
xlabel('X'); ylabel('Y'); zlabel('Z');

h_c3 = subplot(2,2,4); 
scatter3(test_set(:,1), test_set(:,2), test_set(:,3), ...
	3*ones(size(test_set(:,1))), fullcolor_cmixtest, 'filled');
hold on;
scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
	5*ones(n_train, 1), fullcolor_cmixmatch, 'filled');
axis tight;
title('GP Estimated - Variable c');
xlabel('X'); ylabel('Y'); zlabel('Z');

hlinkcone = linkprop([h_quiver1, h_quiver2, h_quiver3], ...
	{'Xlim', 'Ylim', 'Zlim', 'View'});

hlinka = linkprop([h_a3, h_a2, h_a1], {'Xlim', 'Ylim', 'Zlim', 'View'});
hlinkb = linkprop([h_b3, h_b2, h_b1], {'Xlim', 'Ylim', 'Zlim', 'View'});
hlinkc = linkprop([h_c3, h_c2, h_c1], {'Xlim', 'Ylim', 'Zlim', 'View'});
end

%% Results table for path test points

% Calculate estimation at a grid of test points
[a_full_est] = GP_predict(train_set, a_train, [x(:), y(:), z(:)], ...
	cov_funs{1}, fmina) + mean_a;
[b_full_est] = GP_predict(train_set, b_train, [x(:), y(:), z(:)], ...
	cov_funs{1}, fminb) + mean_b;
[c_full_est] = GP_predict(train_set, c_train, [x(:), y(:), z(:)], ...
	cov_funs{1}, fminc) + mean_c;
[a_full_mix] = GP_predict(train_set, a_train, [x(:), y(:), z(:)], ...
	cov_funs{1}, fminab) + mean_a;
[b_full_mix] = GP_predict(train_set, b_train, [x(:), y(:), z(:)], ...
	cov_funs{1}, fminab) + mean_b;
[c_full_mix] = GP_predict(train_set, c_train, [x(:), y(:), z(:)], ...
	cov_funs{1}, fminab) + mean_c;

err_a = sum((a_full_est(:) - u(:)).^2);
err_b = sum((b_full_est(:) - v(:)).^2);
err_c = sum((c_full_est(:) - w(:)).^2);
err_ab = sum((a_full_mix(:) - u(:)).^2) + sum((b_full_mix(:) - v(:)).^2)...
	 + sum((c_full_mix(:) - w(:)).^2);

fprintf('\n----- GP Training: Path Results -----\n');
fprintf(['Method\t\t| %10s\t| %10s\t| %10s\t| %10s\t|| ', ...
    '%s\n'], 'length', 'sigma_f', 'sigma_n', 'nlml', 'Sum Squared Error');
disp(repmat('=', 1, 96));
fprintf('a (alone)\t| %10g\t| %10g\t| %10g\t| %10g\t|| %10g\n', ...
	exp(fmina), nlmla, err_a);
fprintf('b (alone)\t| %10g\t| %10g\t| %10g\t| %10g\t|| %10g\n', ...
	exp(fminb), nlmlb, err_b);
fprintf('c (alone)\t| %10g\t| %10g\t| %10g\t| %10g\t|| %10g\n', ...
	exp(fminc), nlmlc, err_c);
fprintf('Average\t\t| %10g\t| %10g\t| %10g\t| %10g\t|| %10g\n', ...
	(exp(fmina)+exp(fminb)+exp(fminc))/3, nlmla+nlmlb+nlmlc, err_a+err_b+err_c);
disp(repmat('-', 1, 96));
fprintf('Common\t\t| %10g\t| %10g\t| %10g\t| %10g\t|| %10g\n', ...
	exp(fminab), nlmlab, err_ab);

%% Plot comparison full field results
figure(5); clf;
set(gcf, 'name', 'Full field results comparison');
maxV_grid = max(sqrt(u_grid(:).^2 + v_grid(:).^2 + w_grid(:).^2));
maxV_est = max(sqrt(u_est(:).^2 + v_est(:).^2 + w_est(:).^2));
maxV_mix = max(sqrt(u_mix(:).^2 + v_mix(:).^2 + w_mix(:).^2));
maxV = max([maxV_grid, maxV_est, maxV_mix])*1.5;

h_quiver5a = subplot(1,3,1);
daspect([1,1,1]);
h_cone5a = coneplot(x_grid, y_grid, z_grid, u_grid, v_grid, w_grid, ...
	x_grid, y_grid, z_grid, maxV_grid/maxV);
set(h_cone5a, 'EdgeColor', 'none', 'FaceColor', 'blue', 'FaceAlpha', 0.5)
hold on; plot3(x_train, y_train, z_train, 'r+');
axis tight; view(3);
title('Original wind vector field');
xlabel('X'); ylabel('Y');  zlabel('Z');

h_quiver5b = subplot(1,3,2);
daspect([1,1,1]);
h_cone5b = coneplot(x_test, y_test, z_test, u_est, v_est, w_est, ...
	x_test, y_test, z_test, maxV_est/maxV, sqrt(a_std.^2 + b_std.^2 + c_std.^2));
set(h_cone5b, 'EdgeColor', 'none')
hold on; plot3(x_train, y_train, z_train, 'r+');
axis tight; view(3);
title('Seperate HP estimated wind vector field');
xlabel('X'); ylabel('Y');  zlabel('Z');
colorbar('SouthOutside');

h_quiver5c = subplot(1,3,3);
daspect([1,1,1]);
h_cone5c = coneplot(x_test, y_test, z_test, u_mix, v_mix, w_mix, ...
	x_test, y_test, z_test, maxV_mix/maxV, mix_std);
set(h_cone5c, 'EdgeColor', 'none')
hold on; plot3(x_train, y_train, z_train, 'r+');
axis tight; view(3);
title('Common HP estimated wind vector field');
xlabel('X'); ylabel('Y');  zlabel('Z');
colorbar('SouthOutside');

linkprop([h_quiver5a, h_quiver5b, h_quiver5c], ...
	{'Xlim', 'Ylim', 'Zlim', 'View'});

%% Accuracy fields

if exist('x_grid', 'var')
	figure(6); clf;
	set(gcf, 'name', 'Field Accuracy');
	
	h_quiver6a = subplot(1,2,1);
	daspect([1,1,1]);
	h_cone6a = coneplot(x_test, y_test, z_test, u_est, v_est, w_est, ...
		x_test, y_test, z_test, 1, sqrt((u_est - u_grid).^2 + ...
		(v_est - v_grid).^2 + (w_est - w_grid).^2));
	set(h_cone6a, 'EdgeColor', 'none')
	hold on; plot3(x_train, y_train, z_train, 'r+');
	axis tight; view(3);
	title('Seperate HP estimated wind vector field error');
	xlabel('X'); ylabel('Y');  zlabel('Z');
	colorbar('SouthOutside');
	
	h_quiver6b = subplot(1,2,2);
	daspect([1,1,1]);
	h_cone6a = coneplot(x_test, y_test, z_test, u_mix, v_mix, w_mix, ...
		x_test, y_test, z_test, 1, sqrt((u_mix - u_grid).^2 + ...
		(v_mix - v_grid).^2 + (w_mix - w_grid).^2));
	set(h_cone6a, 'EdgeColor', 'none')
	hold on; plot3(x_train, y_train, z_train, 'r+');
	axis tight; view(3);
	title('Common HP estimated wind vector field error');
	xlabel('X'); ylabel('Y');  zlabel('Z');
	colorbar('SouthOutside');
	
	linkprop([h_quiver6a, h_quiver6b], ...
		{'Xlim', 'Ylim', 'Zlim', 'View'});
end

%% SCRAP
% figure(2); clf;
% set(gcf, 'name', 'Seperated hyperparameter solutions');
% h_quiver2 = subplot(2,2,1);
% 
% daspect([1,1,1]);
% h_cone2 = coneplot(x_test, y_test, z_test, u_est, v_est, w_est, ...
% 	x_test, y_test, z_test);
% set(h_cone2, 'EdgeColor', 'none')
% axis tight; hold on;
% title('Estimated wind vector field');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% 
% h_a2 = subplot(2,2,2);
% h_ac2 = coneplot(x_test, y_test, z_test, u_est, zerob, zerob, ...
% 	x_test, y_test, z_test, a_std);
% set(h_ac2, 'EdgeColor', 'none')
% % hold on;
% % scatter3(train_set(:,1), train_set(:,2), train_set(:,3), ...
% % 	5*ones(n_train, 1), fullcolor_amatch, 'filled');
% axis tight;
% title('GP Estimated - Variable a');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% 
% h_b2 = subplot(2,2,3);
% h_bc2 = coneplot(x_test, y_test, z_test, zerob, v_est, zerob, ...
% 	x_test, y_test, z_test, b_std);
% set(h_bc2, 'EdgeColor', 'none')
% axis tight;
% title('GP Estimated - Variable b');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% 
% h_c2 = subplot(2,2,4);
% h_cc2 = coneplot(x_test, y_test, z_test, zerob, zerob, w_est, ...
% 	x_test, y_test, z_test, c_std);
% set(h_cc2, 'EdgeColor', 'none')
% axis tight;
% title('GP Estimated - Variable c');
% xlabel('X'); ylabel('Y'); zlabel('Z');
%
% h_a3 = subplot(2,2,2);
% h_ac3 = coneplot(x_test, y_test, z_test, u_mix, zerob, zerob, ...
% 	x_test, y_test, z_test, mix_std);
% set(h_ac3, 'EdgeColor', 'none')
% axis tight; view(3);
% title('GP Estimated - Variable a');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% 
% h_b3 = subplot(2,2,3);
% h_bc3 = coneplot(x_test, y_test, z_test, zerob, v_mix, zerob, ...
% 	x_test, y_test, z_test, mix_std);
% set(h_bc3, 'EdgeColor', 'none')
% axis tight; view(3);
% title('GP Estimated - Variable b');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% 
% h_c3 = subplot(2,2,4);
% h_cc3 = coneplot(x_test, y_test, z_test, zerob, zerob, w_mix, ...
% 	x_test, y_test, z_test, mix_std);
% set(h_cc3, 'EdgeColor', 'none')
% axis tight; view(3);
% title('GP Estimated - Variable c');
% xlabel('X'); ylabel('Y'); zlabel('Z');