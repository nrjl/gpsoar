%% Setup for running GP over wind data
% clear variables

x_full = -2:0.2:2;
[x_full, y_full] = meshgrid(x_full, x_full);

uv_coord = 1;

%% Flow Selection

% Torus thermal
WW = [0.01, 0.02];
wind_function = @(x, y, t) torus_thermal2d(x-WW(1)*t, y-WW(2)*t+.5);

[a, b] = wind_function(x_full, y_full, 0);

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

% Training points from cross flight
% n_paths = 3;
% n_sections = 2*n_paths+1;
% n_points = 40;
% speed = 0.2;
% 
% x_train = zeros(n_points*n_paths, 1);
% x_onepass = linspace(min(x_full(:)), max(x_full(:)), n_points);
% for i = 1:n_paths
% 	if rem(i, 2)
% 		x_train((i-1)*n_points+1:i*n_points) = x_onepass;
% 	else
% 		x_train((i-1)*n_points+1:i*n_points) = reverse(x_onepass);
% 	end
% end
% 
% deltay = (max(y_full(:)) - min(y_full(:)))/n_sections;
% y_train = zeros(n_points*n_paths, 1);
% for i = 1:n_paths
% 	y_train((i-1)*n_points+1:i*n_points) = ...
% 		linspace((2*i-1)*deltay, (2*i)*deltay, n_points);
% end
% y_train = y_train + min(y_full(:));
% 
% x = [x_train(:), y_train(:)];

y_train = [-2:.05:2]';
x_train = 2*cos(y_train*pi+pi/2);
x_train = [x_train(:);x_train(:)];
y_train = [y_train(:);-y_train(:)];
x = [x_train, y_train];

t = linspace(0, 80, length(x_train));
t = t(:);
[a_obs, b_obs] = wind_function(x_train, y_train, t);

% With noise
a_obs = a_obs + 0.05*randn(size(a_obs));
b_obs = b_obs + 0.05*randn(size(b_obs));

%% Plot training point data
plot(h_quiver1, x(:,1), x(:,2), 'k.');
plot3(h_a1, x(:,1), x(:,2), a_obs, 'k.');
plot3(h_b1, x(:,1), x(:,2), b_obs, 'k.');

%% Do GP analysis

% TEST - Fixed hyperparameters

l_x = 0.3;
l_t = 10;
ab = 2;
sigma_f = 0.2;
sigma_n = 0.1;
cov_funs = {@cressie_ns5, @d_cressie_ns5};
loghyper = ([sigma_f, l_t, l_x, ab, sigma_n]);

cov_funs2 = {@square_expt, @d_square_expt};
loghyper2 = log([l_x, l_t, sigma_f, sigma_n]);

t_star = t(floor(length(t)/2))*ones(size(x_star, 1), 1);

%% GP optimisation

opt = optimset('Display', 'iter', 'GradObj', 'on', 'TolX', 1e-7,...
	'MaxIter', 100);

tta = cputime;
% Separable optimise
[fmin_sep, nlml_sep] = fminunc(@(hyper) ...
	GPt_likelihoodn(x, t, a_obs, b_obs, cov_funs2, hyper), loghyper2, opt);
fprintf(1, '\nSeparable Solution: l_t = %0.5g, lx = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fmin_sep));
[a_sep, a_cov_sep] = ...
	GPt_predict(x, t, a_obs, x_star, t_star, cov_funs2{1}, fmin_sep);
[b_sep, b_cov_sep] = ...
	GPt_predict(x, t, b_obs, x_star, t_star, cov_funs2{1}, fmin_sep);
tta = cputime - tta;

ttb = cputime;
% Non-separable optimise
[fmin_nsep, nlml_nsep] = fminunc(@(hyper) ...
	GPt_likelihoodn(x, t, a_obs, b_obs, cov_funs, hyper), loghyper, opt);
fprintf(1, '\nNon-separable Solution: sf = %0.5g, lt = %0.5g, lx = %0.5g, c = %0.5g, sn = %0.5g\n', ...
	fmin_nsep);
[a_nsep, a_cov_nsep] = ...
	GPt_predict(x, t, a_obs, x_star, t_star, cov_funs{1}, fmin_nsep);
[b_nsep, b_cov_nsep] = ...
	GPt_predict(x, t, b_obs, x_star, t_star, cov_funs{1}, fmin_nsep);
ttb = cputime - ttb;

fprintf('\nEvaluation times: \n %1.4f separable \n %1.4f non-separable\n', tta, ttb);
%% Plot Seperable solution
figure(2); clf;
set(gcf, 'name', 'Seperable covariance function');

% a_covariance = reshape(diag(a_cov_optimum), size(x_full));
% b_covariance = reshape(diag(b_cov_optimum), size(x_full));
a_covariance = reshape(abs(a_cov_sep), size(x_full));
b_covariance = reshape(abs(b_cov_sep), size(x_full));

a_star = reshape(a_sep, size(x_full));
b_star = reshape(b_sep, size(x_full));

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
set(gcf, 'name', 'Non-Seperable covariance function');

% a_covariance = reshape(diag(a_cov_optimum), size(x_full));
% b_covariance = reshape(diag(b_cov_optimum), size(x_full));
a_covariance = reshape(abs(a_cov_nsep), size(x_full));
b_covariance = reshape(abs(b_cov_nsep), size(x_full));

a_star = reshape(a_nsep, size(x_full));
b_star = reshape(b_nsep, size(x_full));

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

h_quiver3 = subplot(2,1,1);
contour(x_full, y_full, a_std);
axis equal; axis tight; hold on;
quiver(x_full, y_full, u_est, v_est, arrowscale);
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');
plot(h_quiver2, x(:,1), x(:,2), 'k.');

h_a3 = subplot(2,2,3);
surf(x_full, y_full, a_star); hold on;
surf(x_full, y_full, a_star + 2*a_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, a_star - 2*a_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimate - {\it u}');
xlabel('X'); ylabel('Y'); zlabel('{\it u}');
plot3(h_a3, x(:,1), x(:,2), a_obs, 'k.');

h_b3 = subplot(2,2,4); 
surf(x_full, y_full, b_star); hold on;
surf(x_full, y_full, b_star + 2*b_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, b_star - 2*b_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimate - {\itv}');
xlabel('X'); ylabel('Y'); zlabel('{\it v}');
plot3(h_b3, x(:,1), x(:,2), b_obs, 'k.');


set([h_a3, h_b3], 'zlim', [-1 1]);
 
 linkaxes([h_quiver1, h_quiver2, h_quiver3]);
 
 hlinka = linkprop([h_a3, h_a2, h_a1], {'Xlim', 'Ylim', 'Zlim', 'View'});
 hlinkb = linkprop([h_b3, h_b2, h_b1], {'Xlim', 'Ylim', 'Zlim', 'View'});


%% List results table

% RMS error
rms_sep = sqrt(mean( (a_sep(:) - a(:)).^2 + (b_sep(:) - b(:)).^2));
rms_nsep = sqrt(mean( (a_nsep(:) - a(:)).^2 + (b_nsep(:) - b(:)).^2));


fprintf('\n----- GP Training results -----\n');
fprintf(['Method\t| l_x\t\t| l_t\t\t|sigma_f\t| sigma_n\t| log(ML)\t|| ', ...
    'RMS Error\n']);
disp(repmat('-', 1, 80));
fprintf('Sep\t\t| %0.5f\t| %0.3f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	exp(fmin_sep),-nlml_sep, rms_sep);
fprintf('NonSep\t| %0.5f\t| %0.3f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	fmin_nsep([3,2,1,5]), -nlml_nsep, rms_nsep);

%% Overplot of all fields
% figure(4); clf; hold on;
% set(4, 'Name', 'Vector field results')
% quiver(x_full, y_full, a, b);
% quiver(x_full, y_full, u_est, v_est, arrowscale);
% quiver(x_full, y_full, u_mix, v_mix, arrowscale_mix);
% legend('Original data', 'Seperate hyperparameter solution', ...
% 	'Common hyperparameter solution');


%% Video
figure(4); clf; set(4, 'position', [500 500 800 320]);
h4 = [0 0 0];
hq4 = [0 0 0];

h4(1) = subplot(1,3,1); hold on;
title('Actual wind field');
[a, b] = wind_function(x_full, y_full, 0);
hq4(1) = quiver(x_full, y_full, a, b);

t_star = 0*ones(size(x_star, 1), 1);

h4(2) = subplot(1,3,2); hold on;
title('Estimate - Seperable covariance function');
[a_sep, a_cov_sep] = GPt_predict(x, t, a_obs, x_star, t_star, cov_funs2{1}, fmin_sep);
[b_sep, b_cov_sep] = GPt_predict(x, t, b_obs, x_star, t_star, cov_funs2{1}, fmin_sep);
a_star = reshape(a_sep, size(x_full));
b_star = reshape(b_sep, size(x_full));
hq4(2) = quiver(x_full, y_full, a_star, b_star);

h4(3) = subplot(1,3,3); hold on;
title('Estimate - Non-seperable covariance function');
[a_nsep, a_cov_nsep] = GPt_predict(x, t, a_obs, x_star, t_star, cov_funs{1}, fmin_nsep);
[b_nsep, b_cov_nsep] = GPt_predict(x, t, b_obs, x_star, t_star, cov_funs{1}, fmin_nsep);
a_star = reshape(a_nsep, size(x_full));
b_star = reshape(b_nsep, size(x_full));
hq4(3) = quiver(x_full, y_full, a_star, b_star);

set(h4, 'xlim', [-2 2], 'ylim', [-2 2]);

hp(1) = plot(h4(1), x(1,1), x(1,2), 'k.');
hp(2) = plot(h4(2), x(1,1), x(1,2), 'k.');
hp(3) = plot(h4(3), x(1,1), x(1,2), 'k.');

t_mov = t;
XX = x;

for ii = 1:length(t_mov)
	delete(hq4); 
	t_star = t_mov(ii)*ones(size(x_star, 1), 1);
	
	[a, b] = wind_function(x_full, y_full, t_mov(ii));
	hq4(1) = quiver(h4(1), x_full, y_full, a, b, 'b');
	
	[a_sep, a_cov_sep] = GPt_predict(x(1:ii,:), t(1:ii), a_obs(1:ii), x_star, t_star, cov_funs2{1}, fmin_sep);
	[b_sep, b_cov_sep] = GPt_predict(x(1:ii,:), t(1:ii), b_obs(1:ii), x_star, t_star, cov_funs2{1}, fmin_sep);
	a_star = reshape(a_sep, size(x_full));
	b_star = reshape(b_sep, size(x_full));
	hq4(2) = quiver(h4(2), x_full, y_full, a_star, b_star, 'b');
	
	a_covariance = reshape(abs(a_cov_sep), size(x_full));
	a_std = sqrt(a_covariance);
	[bah, hc(1)] = contour(h4(2), x_full, y_full, a_std);
	

	[a_nsep, a_cov_nsep] = GPt_predict(x(1:ii,:), t(1:ii), a_obs(1:ii), x_star, t_star, cov_funs{1}, fmin_nsep);
	[b_nsep, b_cov_nsep] = GPt_predict(x(1:ii,:), t(1:ii), b_obs(1:ii), x_star, t_star, cov_funs{1}, fmin_nsep);
	a_star = reshape(a_nsep, size(x_full));
	b_star = reshape(b_nsep, size(x_full));
	hq4(3) = quiver(h4(3), x_full, y_full, a_star, b_star, 'b');
	
	a_covariance = reshape(abs(a_cov_nsep), size(x_full));
	a_std = sqrt(a_covariance);
	[bah, hc(2)] = contour(h4(3), x_full, y_full, -a_std);
		
	set(h4, 'xlim', [-2 2], 'ylim', [-2 2]);
	set(hp, 'xdata', XX(1:ii,1), 'ydata', XX(1:ii,2));
	M(ii) = getframe(gcf);
	
	delete(hc);
end

movie2avi(M, 'H:\Work\Nick\Thesis\code\SeparableCovariance.avi', 'fps', 20, 'Compression','none')