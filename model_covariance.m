% Generate my own spatio-temporal contour plots from wind models
clear variables
close all

thermal_props = [4, 100, 2, 0, 0, -400];
min_point = thermal_props(4:6) - thermal_props(2)*thermal_props(3)*[1,1,2];
max_point = thermal_props(4:6) + thermal_props(2)*thermal_props(3)*[1,1,2];
t_max = 400;

W_steady = [1, 0, -0.5];

W_actual = @(X, t) thermal_field((X - t*W_steady)', thermal_props);

% Random point finder
point_find = @(nn) ones(nn,1)*min_point + rand(nn, 3)*diag((max_point - min_point));

%% Generate covariance data
n_sample = 500;
n_dist = 70; h_vec = linspace(0, 0.5*sqrt(sum((max_point-min_point).^2)), n_dist);
n_time = 70; u_vec = linspace(0, t_max, n_time);

h_scale = thermal_props(2)*thermal_props(3);
u_scale = sqrt(sum(W_steady.^2))*h_scale;

sample_matrix_u = zeros(n_dist, n_time, n_sample);
sample_matrix_v = sample_matrix_u;
sample_matrix_w = sample_matrix_u;
% K = zeros(n_dist, n_time);

for ii = 1:n_dist
	for jj = 1:n_time
		r_points = point_find(n_sample);
		t_points = rand(n_sample, 1)*t_max;
		y_r = W_actual(r_points, t_points);
		
		% Now generate vectors to displaced points (using spherical coords)
		phi_range = rand(n_sample, 1)*2*pi;
		theta_range = rand(n_sample, 1)*pi;
		r_prime = r_points + h_vec(ii)*...
			[sin(theta_range).*cos(phi_range), ...
			 sin(theta_range).*sin(phi_range), ...
			 cos(theta_range)];
		
		t_prime = (round(rand(n_sample, 1))*2 - 1)*u_vec(jj);
		
		y_r_prime = W_actual(r_prime, t_prime);
		
		sample_cov = (y_r.*y_r_prime)';
		
		sample_matrix_u(ii, jj, :) = sample_cov(:,1);
		sample_matrix_v(ii, jj, :) = sample_cov(:,2);
		sample_matrix_w(ii, jj, :) = sample_cov(:,3);
		
	end
end

%% Means
Ku = mean(sample_matrix_u, 3);
Kv = mean(sample_matrix_v, 3);
Kw = mean(sample_matrix_w, 3);

Kus = std(sample_matrix_u, 0, 3);
Kvs = std(sample_matrix_v, 0, 3);
Kws = std(sample_matrix_w, 0, 3);

%%
figure(1); clf;
subplot(2,3,1); imagesc(h_vec, u_vec, Ku); axis tight; %daspect([1 1 1]);
	title('Covariance mean, u'); xlabel('||h||'); ylabel('|u|');
subplot(2,3,2); imagesc(h_vec, u_vec, Kv); axis tight; %daspect([1 1 1]);
	title('Covariance mean, u'); xlabel('||h||'); ylabel('|u|');
subplot(2,3,3); imagesc(h_vec, u_vec, Kw); axis tight; %daspect([1 1 1]);
	title('Covariance mean, v'); xlabel('||h||'); ylabel('|u|');

subplot(2,3,4); imagesc(h_vec, u_vec, Kus); axis tight; %daspect([1 1 1]);
	title('Covariance standard deviation, u'); xlabel('||h||'); ylabel('|u|');
subplot(2,3,5); imagesc(h_vec, u_vec, Kvs); axis tight; %daspect([1 1 1]);
	title('Covariance standard deviation, v'); xlabel('||h||'); ylabel('|u|');
subplot(2,3,6); imagesc(h_vec, u_vec, Kws); axis tight; %daspect([1 1 1]);
	title('Covariance standard deviation, w'); xlabel('||h||'); ylabel('|u|');

%% Crazy plot
figure(3); clf; hold on

for ii = 1:10:n_sample
	mesh(sample_matrix_v(:,:,ii), 'edgealpha', 0, 'facealpha', 0, 'marker', '.')
end

view(3);
		
