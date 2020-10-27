%% Setup for running GP over flight wind data

load Cub_006_data
x_full = x;
y_full = y;

uv_coord = 1;

n_training = 100;

a = u; b = v;

%% Plot original field
figure(1); clf;
quiver(x_full, y_full, u, v);
axis equal; axis tight; hold on;
title('Original wind vector field');
xlabel('X'); ylabel('Y');

% h_a1 = subplot(2,2,3); plot3(x_full, y_full, a, 'b.'); hold on;
% title('a');
% xlabel('X'); ylabel('Y'); zlabel('a');
% 	
% h_b1 = subplot(2,2,4); plot3(x_full, y_full, b, 'b.'); hold on;
% title('b');
% xlabel('X'); ylabel('Y'); zlabel('b');

%% Get training points - Initially a uniform sampling across field

full_set = [x_full(:), y_full(:)];

n_fullset = size(full_set, 1);

X_indexes = floor((rand([1,n_training])*n_fullset));

x = full_set(X_indexes, :);
x_star = full_set;

% Observations
a_obs = a(X_indexes);
b_obs = b(X_indexes);


%% Plot training point data
plot(gcf, h_quiver1, x(:,1), x(:,2), 'r+');

%% Do GP analysis
length_scale = 1;
sigma_f = 1;
sigma_n = 1;
cov_funs = {@square_exp, @d_square_exp};
loghyper = log([length_scale, sigma_f, sigma_n]);

%% GP optimisation

opt = optimset('Display', 'iter', 'GradObj', 'on', 'TolX', 1e-7,...
	'MaxIter', 100);

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

% Group optimise
[fminab, nlmlab] = fminunc(@(hyper) ...
	GP_likelihood2(x, a_obs, b_obs, cov_funs, hyper), loghyper, opt);
fprintf(1, '\nMixed Solution: l = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(fminab(1)), exp(fminab(2)), exp(fminab(3)));
[a_mix, a_cov_mix] = ...
	GP_predict(x, a_obs, x_star, cov_funs{1}, fminab);
[b_mix, b_cov_mix] = ...
	GP_predict(x, b_obs, x_star, cov_funs{1}, fminab);

%% Plot
figure(2); clf;

a_covariance = reshape(diag(a_cov_optimum), size(x_full));
b_covariance = reshape(diag(b_cov_optimum), size(x_full));

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
	max(sqrt(u(:).^2 + v(:).^2));

h_quiver2 = subplot(2,1,1);
contour(x_full, y_full, a_std);
axis equal; axis tight; hold on;
quiver(x_full, y_full, u_est, v_est, arrowscale);
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');
plot(h_quiver2, x(:,1), x(:,2), 'r+');

h_a2 = subplot(2,2,3);
surf(x_full, y_full, a_star); hold on;
surf(x_full, y_full, a_star + 2*a_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, a_star - 2*a_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimated - Variable a');
xlabel('X'); ylabel('Y'); zlabel('a');
plot3(h_a2, x(:,1), x(:,2), a_obs, 'r+');

h_b2 = subplot(2,2,4); 
surf(x_full, y_full, b_star); hold on;
surf(x_full, y_full, b_star + 2*b_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, b_star - 2*b_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimated - Variable b');
xlabel('X'); ylabel('Y'); zlabel('b');
plot3(h_b2, x(:,1), x(:,2), b_obs, 'r+');

%% Plot mixed solution
figure(3); clf;

mix_covariance = reshape(diag(a_cov_mix), size(x_full));

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

arrowscale = max(sqrt(u_mix(:).^2 + v_mix(:).^2))/...
	max(sqrt(u(:).^2 + v(:).^2));

%%
h_quiver3 = subplot(2,1,1);
contour(x_full, y_full, mix_std);
axis equal; axis tight; hold on;
quiver(x_full, y_full, u_est, v_est, arrowscale);
title('Estimated wind vector field');
xlabel('X'); ylabel('Y');
plot(h_quiver3, x(:,1), x(:,2), 'r+');

h_a3 = subplot(2,2,3);
surf(x_full, y_full, a_mix); hold on;
surf(x_full, y_full, a_mix + 2*mix_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, a_mix - 2*mix_std, a_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimated - Variable a');
xlabel('X'); ylabel('Y'); zlabel('a');
plot3(h_a3, x(:,1), x(:,2), a_obs, 'r+');

h_b3 = subplot(2,2,4); 
surf(x_full, y_full, b_mix); hold on;
surf(x_full, y_full, b_mix + 2*mix_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(x_full, y_full, b_mix - 2*mix_std, b_colourcorrect, ...
	'FaceAlpha', 0.5, 'EdgeColor', 'none');
title('GP Estimated - Variable b');
xlabel('X'); ylabel('Y'); zlabel('b');
plot3(h_b3, x(:,1), x(:,2), b_obs, 'r+');

%%

linkaxes([h_quiver1, h_quiver2, h_quiver3]);

hlinka = linkprop([h_a3, h_a2, h_a1], {'Xlim', 'Ylim', 'Zlim', 'View'});
hlinkb = linkprop([h_b3, h_b2, h_b1], {'Xlim', 'Ylim', 'Zlim', 'View'});


%% List results table
err_a = sum((a_optimum(:) - u(:)).^2);
err_b = sum((b_optimum(:) - v(:)).^2);
err_ab = sum((a_mix(:) - u(:)).^2) + sum((b_mix(:) - v(:)).^2);

fprintf('\n----- GP Training results -----\n');
fprintf(['Method\t\t| length\t| sigma_f\t| sigma_n\t| nlml\t\t|| ', ...
    'Sum Squared Error\n']);
fprintf('a (alone)\t| %0.5g\t| %0.5g\t| %0.5g\t| %0.5g\t|| %0.5g\n', ...
	exp(fmina), nlmla, err_a);
fprintf('b (alone)\t| %0.5g\t| %0.5g\t| %0.5g\t| %0.5g\t|| %0.5g\n', ...
	exp(fminb), nlmlb, err_b);
fprintf('Grouped\t\t| %0.5g\t| %0.5g\t| %0.5g\t| %0.5g\t|| %0.5g\n', ...
	exp(fminab), nlmlab, err_ab);