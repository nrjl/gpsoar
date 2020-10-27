%% Setup for running GP over wind data
% clear variables

x_full = -2:0.2:2;
[x_full, y_full] = meshgrid(x_full, x_full);

uv_coord = 1;

%% Flow Selection

% % Inviscid ideal flow field
% flow_c = {@vortex, [0, 1, 1]; 
% 		  @vortex, [0, -1, -1];
% 		  @source, [1, 0, 1];
% 		  @stream, [1, pi]};
% 
% wind_function = @(x, y) make_field(flow_c, x, y);

% Torus thermal
wind_function = @torus_thermal2d;

[a, b] = wind_function(x_full, y_full);

%% Plot original field
figure(1); clf;
set(gcf, 'name', 'Original system');
h_quiver1 = subplot(2,1,1); quiver(x_full, y_full, a, b);
axis equal; axis tight; hold on;
title('Original wind vector field');
xlabel('X'); ylabel('Y');

h_a1 = subplot(2,2,3); surf(x_full, y_full, a); hold on;
title('\it u');
xlabel('X'); ylabel('Y'); zlabel('\it u');
	
h_b1 = subplot(2,2,4); surf(x_full, y_full, b); hold on;
title('\it v');
xlabel('X'); ylabel('Y'); zlabel('\it v');

full_set = [x_full(:), y_full(:)];
x_star = full_set;

n_fullset = size(full_set, 1);
set([h_a1, h_b1], 'zlim', [-1 1]);

%% Get training points 
% % Initially a uniform sampling across field
% % spacing = floor(n_fullset/n_training);
% % X_indexes = floor(spacing/2):spacing:n_fullset;
% n_training = 30;
% X_indexes = floor((rand([1,n_training])*n_fullset));
% 
% x = full_set(X_indexes, :);
% 
% % Observations (no noise)
% a_obs = a(X_indexes)';
% b_obs = b(X_indexes)';

% Training points from cross flight
n_paths = 3;
n_sections = 2*n_paths+1;
n_points = 40;

x_train = zeros(n_points*n_paths, 1);
x_onepass = linspace(min(x_full(:)), max(x_full(:)), n_points);
for i = 1:n_paths
	if rem(i, 2)
		x_train((i-1)*n_points+1:i*n_points) = x_onepass;
	else
		x_train((i-1)*n_points+1:i*n_points) = reverse(x_onepass);
	end
end

deltay = (max(y_full(:)) - min(y_full(:)))/n_sections;
y_train = zeros(n_points*n_paths, 1);
for i = 1:n_paths
	y_train((i-1)*n_points+1:i*n_points) = ...
		linspace((2*i-1)*deltay, (2*i)*deltay, n_points);
end
y_train = y_train + min(y_full(:));

[a_obs, b_obs] = wind_function(x_train, y_train);
x = [x_train(:), y_train(:)];

% With noise
a_obs = a_obs + 0.05*randn(size(a_obs));
b_obs = b_obs + 0.05*randn(size(b_obs));

%% Plot training point data
plot(h_quiver1, x(:,1), x(:,2), 'k.');
plot3(h_a1, x(:,1), x(:,2), a_obs, 'k.');
plot3(h_b1, x(:,1), x(:,2), b_obs, 'k.');

%% Do GP analysis

% TEST - Fixed hyperparameters

length_scale = 1;
sigma_f = 0.1;
sigma_n = 0.1;
cov_funs = {@square_exp, @d_square_exp};
loghyper = log([length_scale, sigma_f, sigma_n]);

%% GP optimisation


opt = optimset('Display', 'iter', 'GradObj', 'on', 'TolX', 1e-7,...
	'MaxIter', 100);

tta = cputime;
% a optimise
[fmina, nlmla] = fminunc(...
	@(hyper) GP_likelihood(x, a_obs, cov_funs, hyper), loghyper, opt);
fprintf(1, '\nSolution a: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fmina(1)), exp(fmina(2)), exp(fmina(3)));
[a_optimum, a_cov_optimum] = ...
	GP_predict(x, a_obs, x_star, cov_funs{1}, fmina);

% b optimise
[fminb, nlmlb] = fminunc(...
	@(hyper) GP_likelihood(x, b_obs, cov_funs, hyper), loghyper, opt);

fprintf(1, '\nSolution b: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fminb(1)), exp(fminb(2)), exp(fminb(3)));
[b_optimum, b_cov_optimum] = ...
	GP_predict(x, b_obs, x_star, cov_funs{1}, fminb);
tta = cputime - tta;

ttb = cputime;
% Group optimise
[fminab, nlmlab] = fminunc(@(hyper) ...
	GP_likelihoodn(x, a_obs, b_obs, cov_funs, hyper), loghyper, opt);
fprintf(1, '\nMixed Solution: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fminab(1)), exp(fminab(2)), exp(fminab(3)));
[a_mix, a_cov_mix] = ...
	GP_predict(x, a_obs, x_star, cov_funs{1}, fminab);
[b_mix, b_cov_mix] = ...
	GP_predict(x, b_obs, x_star, cov_funs{1}, fminab);
ttb = cputime - ttb;

fprintf('\nEvaluation times: \n %1.4f separate \n %1.4f shared\n', tta, ttb);
%% Plot
figure(2); clf;
set(gcf, 'name', 'Seperated hyperparameter solutions');

% a_covariance = reshape(diag(a_cov_optimum), size(x_full));
% b_covariance = reshape(diag(b_cov_optimum), size(x_full));
a_covariance = reshape(a_cov_optimum, size(x_full));
b_covariance = reshape(b_cov_optimum, size(x_full));

a_star = reshape(a_optimum, size(x_full));
b_star = reshape(b_optimum, size(x_full));

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
	max(sqrt(a(:).^2 + b(:).^2));

h_quiver2 = subplot(2,1,1);
contour(x_full, y_full, a_std);
axis equal; axis tight; hold on;
quiver(x_full, y_full, u_est, v_est, arrowscale);
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');
plot(h_quiver2, x(:,1), x(:,2), 'k.');

h_a2 = subplot(2,2,3);
surf(x_full, y_full, a_star); hold on;
surf(x_full, y_full, a_star + 2*a_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, a_star - 2*a_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimate - {\it u}');
xlabel('X'); ylabel('Y'); zlabel('{\it u}');
plot3(h_a2, x(:,1), x(:,2), a_obs, 'k.');

h_b2 = subplot(2,2,4); 
surf(x_full, y_full, b_star); hold on;
surf(x_full, y_full, b_star + 2*b_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, b_star - 2*b_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimate - {\itv}');
xlabel('X'); ylabel('Y'); zlabel('{\it v}');
plot3(h_b2, x(:,1), x(:,2), b_obs, 'k.');


set([h_a2, h_b2], 'zlim', [-1 1]);

%% Plot mixed solution
figure(3); clf;
set(gcf, 'name', 'Common hyperparameter solution');

% mix_covariance = reshape(diag(a_cov_mix), size(x_full));
mix_covariance = reshape(a_cov_mix, size(x_full));

a_mix = reshape(a_mix, size(x_full));
b_mix = reshape(b_mix, size(x_full));

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

arrowscale_mix = max(sqrt(u_mix(:).^2 + v_mix(:).^2))/...
	max(sqrt(a(:).^2 + b(:).^2));

%%
h_quiver3 = subplot(2,1,1);
contour(x_full, y_full, mix_std);
axis equal; axis tight; hold on;
quiver(x_full, y_full, u_mix, v_mix, arrowscale_mix);
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');
plot(h_quiver3, x(:,1), x(:,2), 'k.');

h_a3 = subplot(2,2,3);
surf(x_full, y_full, a_mix); hold on;
surf(x_full, y_full, a_mix + 2*mix_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, a_mix - 2*mix_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimate - {\it u}');
xlabel('X'); ylabel('Y'); zlabel('{\it u}');
plot3(h_a3, x(:,1), x(:,2), a_obs, 'k.');

h_b3 = subplot(2,2,4); 
surf(x_full, y_full, b_mix); hold on;
surf(x_full, y_full, b_mix + 2*mix_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, b_mix - 2*mix_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimate - {\it v}');
xlabel('X'); ylabel('Y'); zlabel('{\it u}');
plot3(h_b3, x(:,1), x(:,2), b_obs, 'k.');


set([h_a3, h_b3], 'zlim', [-1 1]);

%%

linkaxes([h_quiver1, h_quiver2, h_quiver3]);

hlinka = linkprop([h_a3, h_a2, h_a1], {'Xlim', 'Ylim', 'Zlim', 'View'});
hlinkb = linkprop([h_b3, h_b2, h_b1], {'Xlim', 'Ylim', 'Zlim', 'View'});


%% List results table
err_a = sum((a_optimum(:) - a(:)).^2);
err_b = sum((b_optimum(:) - b(:)).^2);
err_ab = sum((a_mix(:) - a(:)).^2) + sum((b_mix(:) - b(:)).^2);

fprintf('\n----- GP Training results -----\n');
fprintf(['Method\t\t| length\t| sigma_f\t| sigma_n\t| nlml\t\t|| ', ...
    'Sum Squared Error\n']);
disp(repmat('-', 1, 80));
fprintf('u (alone)\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	exp(fmina), nlmla, err_a);
fprintf('v (alone)\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	exp(fminb), nlmlb, err_b);
fprintf('Grouped\t\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	exp(fminab), nlmlab, err_ab);

% RMS error
rms_a = sqrt(mean((a_optimum(:) - a(:)).^2));
rms_b = sqrt(mean((b_optimum(:) - b(:)).^2));
rms_aob = sqrt(mean( (a_optimum(:) - a(:)).^2 + (b_optimum(:) - b(:)).^2));

rms_ab = sqrt(mean( (a_mix(:) - a(:)).^2 + (b_mix(:) - b(:)).^2));



fprintf('\n----- GP Training results -----\n');
fprintf(['Method\t\t| length\t| sigma_f\t| sigma_n\t| log(ML)\t|| ', ...
    'RMS Error\n']);
disp(repmat('-', 1, 80));
fprintf('u (alone)\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	exp(fmina),-nlmla, rms_a);
fprintf('v (alone)\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	exp(fminb), -nlmlb, rms_b);
fprintf('uv (total)\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	[0,0,0], -nlmla - nlmlb, rms_aob);
fprintf('Grouped\t\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	exp(fminab), -nlmlab, rms_ab);

%% Overplot of all fields
figure(4); clf; hold on;
set(4, 'Name', 'Vector field results')
quiver(x_full, y_full, a, b);
quiver(x_full, y_full, u_est, v_est, arrowscale);
quiver(x_full, y_full, u_mix, v_mix, arrowscale_mix);
legend('Original data', 'Seperate hyperparameter solution', ...
	'Common hyperparameter solution');