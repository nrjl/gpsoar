% Spatio-temporal ordinary kriging

% Sensor noise
sigma_s = 0.05;

% Generate wind function, must be a function of time and space
thermal_props = [4, 100, 3, 0, 0, -200];
W_steady = [2, 0, -0.5];
W_actual = @(X, t) torus_thermal3((X - t'*W_steady)', thermal_props);

% Training points

% -- Helical path -- %
ri = 75;
ro = 150;
zb = -100; zt = -300;
n_cycles = 4;
n_samples = 500;

[x_train, y_train, z_train] = ...
 	helical_path(ri, ro, zb, zt, n_cycles, n_samples);
X_train = [x_train(:), y_train(:), z_train(:)];
t_train = linspace(0, 50, numel(x_train));

W_train = W_actual(X_train, t_train) + sigma_s*randn(size(X_train));

% Choose (semi-)variogram
vario = @vario1;

% Optimise hyper-parameters (using weighted least squares WLS)
[hyper_min, wls_val] = fminunc(...
	@(hyper) wls(X_train, t_train, W_train, vario, hyper), hyper, opt);

% Show solution
