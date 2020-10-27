%% TSHIFT SIMPLE EXAMPLE

% 2D fields

% Static thermal
% thermal_props = [1, 1, 3, -0.5, 0, -2];
thermal_props = [1, 1, 3, -0.5, 0, 0];
W_steady = [0.00, 0.00];
core_start = 0.1;
core_rate = 0.00;

W_actual = @(X, t) thermal_props(1)*(core_start*eye(numel(t)) + core_rate*diag(t))*torus_thermal2d(X(:,1)-W_steady(1)*t, X(:,2)-W_steady(2)*t, ...
	[1, thermal_props(2:end)]); %+repmat(W_steady, [numel(t), 1]);

t_final = 100;

core_max = core_start+core_rate*t_final;
core_scale = @(current_t) (core_start + core_rate*current_t)/core_max;

future_data = false;


%% Plot field at t = 0
xx = linspace(-2, 2, 15)';
yy = linspace(-2, 2, 15)';
[xx, yy] = meshgrid(xx, yy);

WW = W_actual([xx(:), yy(:)], zeros(numel(xx), 1));
UU = reshape(WW(:,1), size(xx));
VV = reshape(WW(:,2), size(yy));
figure(1); clf;
quiver(xx, yy, UU, VV);


%% Observations
% % Transection paths
% n_paths = 5;
% n_sections = 2*n_paths+1;
% n_points = 30;
% 
% x_train = zeros(n_points*n_paths, 1);
% x_onepass = linspace(min(xx(:)), max(xx(:)), n_points);
% for i = 1:n_paths
% 	if rem(i, 2)
% 		x_train((i-1)*n_points+1:i*n_points) = x_onepass;
% 	else
% 		x_train((i-1)*n_points+1:i*n_points) = reverse(x_onepass);
% 	end
% end
% 
% deltay = (max(yy(:)) - min(yy(:)))/n_sections;
% y_train = zeros(n_points*n_paths, 1);
% for i = 1:n_paths
% 	y_train((i-1)*n_points+1:i*n_points) = ...
% 		linspace((2*i-1)*deltay, (2*i)*deltay, n_points);
% end
% y_train = y_train + min(yy(:));
% X = [x_train(:), y_train(:)];

% Circular path
n_circles = 3;
crad = 1;
n_points = 150;
ctheta = linspace(0, n_circles*2*pi, n_points);

X = [crad*cos(ctheta(:)), crad*sin(ctheta(:))];
t = linspace(0, t_final, size(X, 1))';
Y = W_actual(X, t);

% Add noise
Y = Y+0.01*randn(size(Y));


%% Field video
% figure(1); clf;
% h_q = quiver(xx, yy, UU, VV, core_scale(t(1))); axis tight; hold on;
% h_p = plot(X(1,1), X(1,2), 'rx');
% 
% for ii = 2:length(t)
% 	WW = W_actual([xx(:), yy(:)], t(ii)*ones(numel(xx), 1));
% 	UU = reshape(WW(:,1), size(xx));
% 	VV = reshape(WW(:,2), size(yy));
% 	
% 	set(h_q, 'UData', UU, 'VData', VV, 'AutoScaleFactor', core_scale(t(ii)));
% 	set(h_p, 'XData', X(1:ii, 1), 'YData', X(1:ii, 2));
% 	M_field(ii) = getframe;
% end



%% Plot distance metrics
% Plot observations
figure(2); clf;
% scaleyscale = sqrt(max(sum(Y.^2, 2)))/sqrt(max(UU(:).^2 + VV(:).^2));
subplot(2,2,1);
quiver(X(:,1), X(:,2), Y(:,1), Y(:,2), 1, 'r');

d_xx = square_dist(X, X);
d_tt = square_dist(t, t);

X_prime = X + diag(t(1)-t,0)*Y;
% d_pp = projected_point_distance(X, t, Y, X, t(end));
d_pp = square_dist(X_prime, X_prime);

subplot(2,2,2);
imagesc(d_xx);

subplot(2,2,3);
imagesc(d_tt);

subplot(2,2,4);
imagesc(d_pp);

%% Plot Covariances
% Square-exponential spatial covariance
length_scale = 1;	% Spatial units
plength_scale = 20;	% Projected points length scale
sigma_f = 1;
sigma_n = 0.02;		% NOTE: THIS COVERS TIME NOISE AS WELL!! Does it??

% Exponential temporal covariance
lag_scale = 50;		% Temporal units

opt = optimset('Display', 'iter', 'GradObj', 'off', 'TolX', 1e-5,...
	'MaxIter', 100);

[pfmin, pnlml] = fminunc(@(hyper) GPtp_likelihood(X, t, Y, ...
	{@projected_square_expt}, hyper), ...
	log([length_scale, lag_scale, plength_scale, sigma_f, sigma_n]), opt);
length_scale = exp(pfmin(1)); lag_scale = exp(pfmin(2));
plength_scale = exp(pfmin(3)); sigma_f = exp(pfmin(4));
sigma_n = exp(pfmin(5));

% Display results
disp('----- Results: Projected-point square exponential -----');
fprintf(1, '\nl = %0.5g, lt = %0.5g, lp = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(pfmin));

Kxx = square_exp(X, X, log([length_scale, sigma_f]), d_xx);
Kxt = square_expt(X, t, X, t, log([length_scale, lag_scale, sigma_f]), ...
	d_xx, d_tt);
Kxp = projected_square_expt(X, t, Y, X, t, ...
	log([length_scale, lag_scale, plength_scale, sigma_f]), d_xx, d_tt, d_pp);

figure(3); clf;
subplot(2,2,1);
imagesc(Kxx);

subplot(2,2,2);
imagesc(Kxt);

subplot(2,2,3);
imagesc(Kxp);

%% sq_expt solution
cov_funsa = {@square_expt, @d_square_expt};
loghypera = log([length_scale, lag_scale, sigma_f, sigma_n]);

% Gaussian mean function
% NOTE: B values are standard deviation, which are squared to variance to
% prevent negative values.
mean_fun = @(X) ones([1, size(X, 1)]);
order_h = 1;

supercov_funs = {@square_expt};
W_mean = mean(Y);
superloghyper = [loghypera, ([W_mean(1), 1, W_mean(2), 1])];
supern_hyper = [4, 4];

[gfmin, gnlml] = fminunc(@(hyper) ...
	GPt_likelihoodn_gaussmean(X, t, Y(:,1), Y(:,2), ...
	supercov_funs, mean_fun, hyper, supern_hyper), superloghyper, opt);
% gfmin = positivify(gfmin, supern_hyper, order_h);
[K_inv2] = GPt_predict_gaussmean(X, t, Y(:,1), [], [], supercov_funs{1}, mean_fun, gfmin, supern_hyper);


% Display results
disp('----- Results: Spatio-temporal squared exponential and gaussian mean function -----');
fprintf(1, '\nsq_expt function Solution: l = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', ...
	exp(gfmin(1)), exp(gfmin(2)), exp(gfmin(3)), exp(gfmin(4)));
fprintf(1, 'Wx = %0.5g, sigma_x = %0.5g, Wy = %0.5g, sigma_y = %0.5g\n',...
	gfmin(5:end));

%% Covariance movie (craziness)
figure(4); clf; set(4, 'position', [100 100 960 960]);
Ktt = square_exp(t, t, log([lag_scale, sigma_f]), d_tt);

if future_data
	data_step = size(X,1);
else
	data_step = 1;
end

h4a = subplot(2,2,1);hold on;
h_q = quiver(h4a, xx, yy, UU, VV,  core_scale(t(1))); 
set(h_q, 'autoscale', 'off');
h_p = plot(h4a, X(1,1), X(1,2), 'rx');
h_p2 = scatter(h4a, X_prime(1:data_step,1), X_prime(1:data_step,2), [], Ktt(1:data_step,1));
% h_p2 = plot(h4a, X_prime(:,1), X_prime(:,2), 'gx');
title('Field and estimate points');


% Evaluate prediction
[W_est, cov_est] = GPtp_predict(X(1:data_step,:), t(1:data_step), Y(1:data_step,:), [xx(:), yy(:)], 0, log([length_scale, lag_scale, plength_scale, sigma_f]));
U_est = reshape(W_est(:,1), size(xx));
V_est = reshape(W_est(:,2), size(yy));
h_q_est = quiver(h4a, xx, yy, U_est, V_est, core_scale(t(1)), 'g');
set(h_q_est, 'autoscale', 'off'); axis tight; 
[HC, h_contour] = contour(xx, yy, reshape(cov_est, size(xx)));

h4b = subplot(2,2,2);
hdpp = imagesc(d_pp);
title('Projected point distance matrix');

h4c = subplot(2,2,3);
hspp = imagesc(square_exp(X_prime, X_prime, log([plength_scale, sigma_f]), d_pp));
title('Square exponential covariance based on projected distances');

h4d = subplot(2,2,4);
hsxp = imagesc(projected_square_expt(X, t, Y, X, t, ...
	log([length_scale, lag_scale, plength_scale, sigma_f]), d_xx, d_tt, d_pp));
title('Full spatio-termporal sq exp covariance based on d_{xx}, d_{tt} and d_{pp}');

M2(1) = getframe(4);

figure(5); clf; 
h_error = axes; set(h_error, 'Nextplot', 'add'); set(h_error, 'XLim', [0 length(t)]);

for ii = 2:length(t)
	if future_data
		data_step = length(t);
	else
		data_step = ii;
	end
	
	WW = W_actual([xx(:), yy(:)], t(ii)*ones(numel(xx), 1));
	UU = reshape(WW(:,1), size(xx));
	VV = reshape(WW(:,2), size(yy));
	
	X_prime = X(1:data_step,:) + diag(t(ii)-t(1:data_step),0)*Y(1:data_step,:);
	d_pp2 = square_dist(X_prime, X_prime);
	Kpp2 = square_exp(X_prime, X_prime, log([plength_scale, sigma_f]), d_pp2);
	Kxp2 = projected_square_expt(X(1:data_step,:), t(1:data_step), Y(1:data_step,:), X(1:data_step,:), t(1:data_step), ...
		pfmin, d_xx(1:data_step,1:data_step), d_tt(1:data_step,1:data_step), d_pp2(1:data_step,1:data_step));
	
	K_mat = Kxp2 + sigma_n^2*eye(size(X_prime,1)); % Add noise
	L = chol(K_mat, 'lower');
	K_inv = (L'\(L\eye(size(L))));	

	% Evaluate prediction
	[W_est, cov_est] = GPtp_predict(X(1:data_step,:), t(1:data_step), Y(1:data_step,:), [xx(:), yy(:)], t(ii), ...
		log([length_scale, lag_scale, plength_scale, sigma_f]), K_inv);
	U_est = reshape(W_est(:,1), size(xx));
	V_est = reshape(W_est(:,2), size(yy));
	
	
	set(h_q, 'UData', UU, 'VData', VV);
	set(h_q_est, 'UData', U_est, 'VData', V_est);
	set(h_contour, 'ZData', reshape(cov_est, size(xx)));
	
	set(h_p, 'XData', X(1:ii, 1), 'YData', X(1:ii, 2));
	set(h_p2, 'XData', X_prime(1:data_step, 1), 'YData', X_prime(1:data_step, 2), 'CData', Ktt(ii,1:data_step));
	
 	set(hdpp, 'CData', d_pp2);
	set(hspp, 'CData', Kpp2);
	set(hsxp, 'Cdata', Kxp2);
 	
	M2(ii) = getframe(4);
	
	plot(h_error, ii, mean((UU(:) - U_est(:)).^2+(VV(:) - V_est(:)).^2), 'b.');
		
	[W_est2, V_est2] = GPt_predictn_gaussmean(X(1:data_step,:), t(1:data_step),  Y(1:data_step,:), ...
		[xx(:), yy(:)], t(ii)*ones(numel(xx(:)), 1), supercov_funs{1}, mean_fun, gfmin, supern_hyper);
	plot(h_error, ii, mean((UU(:) - W_est2(:,1)).^2 + (VV(:) - W_est2(:,2)).^2), 'rx');	
end

legend(h_error, {'Projected point MSE', 'Square exponential MSE'});

%% Make movie
movie2avi(M2, 'movies\GPtp_movie.avi', 'compression', 'none', 'fps', 10)