%% Setup for running GP over flight wind data
clear all
close all

%%
% Get training data from log file. Log contains x and y coordinates stored
% as x and y vectors and corresponding u and v wind components
load Cub_006_data
dstart = 1;
dstop = numel(x);
x = x(dstart:dstop);
y = y(dstart:dstop);
u = u(dstart:dstop);
v = v(dstart:dstop);


train_skip = 30;		% Skip number of training set (ie use every nth data point)

x_train = x(1:train_skip:end);
y_train = y(1:train_skip:end);
n_train = numel(x_train);

train_set = [x_train, y_train];

a_train = u(1:train_skip:end);
b_train = v(1:train_skip:end);

% Uniform test field
nx_test = 20;
ny_test = 15;

uv_coord = 1;


%% Plot original field (which is effectively training data)
time_color = jet(n_train);

figure(1); clf;
set(gcf, 'Name', 'Original data');
h_quiver1 = subplot(2,1,1);
quiver(x, y, u, v);
axis equal; axis tight; hold on;
title('Original wind vector field');
xlabel('X'); ylabel('Y');

h_a1 = subplot(2,2,3); hold on;
% plot3(train_set(:,1), train_set(:,2), a_train, 'r+');
scatter3(train_set(:,1), train_set(:,2), a_train, 4*ones(n_train, 1), time_color, 'filled');
axis tight; view(3);
title('a');
xlabel('X'); ylabel('Y'); zlabel('a');
	
h_b1 = subplot(2,2,4); hold on;
% plot3(train_set(:,1), train_set(:,2), b_train, 'r+');
scatter3(train_set(:,1), train_set(:,2), b_train, 4*ones(n_train, 1), time_color, 'filled');
axis tight; view(3);
title('b');
xlabel('X'); ylabel('Y'); zlabel('b');

drawnow;
% return;

%% Get test points - Uniform sampling across field

x_test = linspace(min(x_train), max(x_train), nx_test);
y_test = linspace(min(y_train), max(y_train), ny_test);
[x_test, y_test] = meshgrid(x_test, y_test);

test_set = [x_test(:), y_test(:)];

n_test_set = size(test_set, 1);


%% Do GP analysis
length_scale = 1000;
sigma_f = 1000;
sigma_n = 1;
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

% b optimise
[fminb, nlmlb] = fminunc(...
	@(hyper) GP_likelihood(train_set, b_train, cov_funs, hyper), loghyper, opt);

fprintf(1, '\nSolution b: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fminb(1)), exp(fminb(2)), exp(fminb(3)));
[b_optimum, b_cov_optimum] = ...
	GP_predict(train_set, b_train, test_set, cov_funs{1}, fminb);

% Group optimise
[fminab, nlmlab] = fminunc(@(hyper) ...
	GP_likelihood2(train_set, a_train, b_train, cov_funs, hyper), loghyper, opt);
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

%%
figure(2); clf;
set(gcf, 'Name', 'Separate solutions');
h_quiver2 = subplot(2,1,1);
contour(x_test, y_test, a_std);

axis equal; axis tight; hold on;
quiver(x_test, y_test, u_est, v_est, arrowscale);
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');
plot(h_quiver2, train_set(:,1), train_set(:,2), 'r+');

h_a2 = subplot(2,2,3);
surf(x_test, y_test, a_star); hold on;
surf(x_test, y_test, a_star + 2*a_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_test, y_test, a_star - 2*a_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimated - Variable a');
xlabel('X'); ylabel('Y'); zlabel('a');
plot3(h_a2, train_set(:,1), train_set(:,2), a_train, 'r+');

h_b2 = subplot(2,2,4); 
surf(x_test, y_test, b_star); hold on;
surf(x_test, y_test, b_star + 2*b_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_test, y_test, b_star - 2*b_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimated - Variable b');
xlabel('X'); ylabel('Y'); zlabel('b');
plot3(h_b2, train_set(:,1), train_set(:,2), b_train, 'r+');

%% Plot mixed solution
figure(3); clf;

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

%%
h_quiver3 = subplot(2,1,1);
contour(x_test, y_test, mix_std);
axis equal; axis tight; hold on;
quiver(x_test, y_test, u_est, v_est, arrowscale);
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');
plot(h_quiver3, train_set(:,1), train_set(:,2), 'r+');

h_a3 = subplot(2,2,3);
surf(x_test, y_test, a_mix); hold on;
surf(x_test, y_test, a_mix + 2*mix_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_test, y_test, a_mix - 2*mix_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimated - Variable a');
xlabel('X'); ylabel('Y'); zlabel('a');
plot3(h_a3, train_set(:,1), train_set(:,2), a_train, 'r+');
axis tight;

h_b3 = subplot(2,2,4); 
surf(x_test, y_test, b_mix); hold on;
surf(x_test, y_test, b_mix + 2*mix_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_test, y_test, b_mix - 2*mix_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimated - Variable b');
xlabel('X'); ylabel('Y'); zlabel('b');
plot3(h_b3, train_set(:,1), train_set(:,2), b_train, 'r+');
axis tight;

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