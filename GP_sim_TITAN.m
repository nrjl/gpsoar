% Use GPs to predict wind, then use RRT-like planner to plan energy gain
% path 
clearvars -except full_turb
addpath static_functions wind_functions

%% SETUP SIMULATION
global Cd0 S AR e m Nmax Nmin CL_max dphi_dt_max GR_approx V_stall phi_max gamma_max
% initialise_SBXC 
global g pitch_coeff
global x_limits y_limits z_limits
global K_explore
global zt1999

load yelle1999.mat zt1999

g = 1.35;

GPt = false;
retrain = false;

% --- Aircraft parameters --- %
Cd0		= 0.02;		% Parasitic drag coefficient
S		= 4.0;		% Wing reference area
b		= 6.25;			% Wing span
AR		= b*b/S;		% Wing aspect ratio
e		= 0.9;			% Oswald's efficiency factor
m		= 100;			% Vehicle mass
Nmax	= 1.7;			% Maximum load factor (positive)
Nmin	= 0.3;			% Minimum load factor (negative)
CL_max	= 1.0;			% Maximum lift coefficient
dphi_dt_max = 30*pi/180;% Maximum roll rate (rad/s)
phi_max = 45*pi/180;	% Max bank
gamma_max = 40*pi/180;	% Max 

max_climb = 50*pi/180;	% Maximum air relative climb angle (rad)
GR_approx = 12;
V_stall = 6;			% Stall speed (m/s)
u_20 = 0;

K_explore = 0.5;

% --- Figure properties --- %
fig_height = 480;
fig_width = 640;
movie_on = false;

% --- Wind field --- %
thermals = [1, 100, 4, 150, 50, -200; 
			1, 100, 4, -150, -50, -200];
shear_props = [-200, -170, 3, 0, 1, 3;
			   -230, -200, -3, 0, 1, 3];

wave_props = [-2, 1000];

t1 = @(X, P) thermal_field(X, P);
s1 = @(X, P) cos_profile3(X, P);		% dWx_dz = [1,3]
w1 = @(X, P) sinusoidal_wind3(X, P);

windspeed = [0, 0, 0]';

% W_actual = @(X) multi_field3(X, s1, shear_props(2,:));
% W_actual = @(X) multi_field3(X, t1, thermals, s1, shear_props(1,:));
W_actual = @(X, t) multi_field3(X, w1, wave_props);

wind_noise = 0.05;			% std of noise on wind observations

max_wind = 3;

%% SETUP GP PROPERTIES

if GPt == false		% Spatial only
	loghyper0 = log([270, 0.8, 0.05]);
	sigma_hyper = log([100, 0.5, 0.5]);
	% loghyper0 = [4.3049   -0.2166   -2.9661];
	cov_funs = {@square_exp, @d_square_exp, @dsquare_exp_dx};
	
else				% Spatio-temporal
% 	loghyper0 = log([30, 100, 0.5, 0.09]);
% 	sigma_hyper = log([15, 30, 0.5, 0.5]);
% 	cov_funs = {@square_expt, @d_square_expt, @dsquare_exp_dxt};
	
	loghyper0 = [log([50, 80, 80, 1]), 0, 0, 0, 0,log(0.5)];
	sigma_hyper = log([15, 20, 20, 1, 0.5, 0.5, 0.5, 0.5, 0.5]);
	cov_funs = {@square_exptdrift, [], @dsquare_exptdrift_dxt};
end

loghyper = loghyper0;

save_hypers = zeros(100, size(loghyper0,2)+1);
save_hypers(1,:) = [0, loghyper0];
hyper_count = 2;

% Observation set
max_obs = 150;
freq_obs = 1;		% Observation frequency

% Maximum target box
x_limits = [0, 1500];
y_limits = [-1000, 1000];
z_limits = [-25000, -24800];

%% INITIAL FLIGHT (MANUAL OR AN A-PRIORI DATA SET)
% Intial wind estimate from manual control
p0 = linspace(0, -pi, 7);
x0 = [linspace(400, 0, 20), 50*sin(p0(2:end))];
y0 = [linspace(-100, 160, 20), 110+50*cos(p0(2:end))]; 
z0 = linspace(-24920, -24900, numel(x0));

% x0 =[-10, 0]; y0 =[0, 0]; z0 = [-205 -200];

X_init = [x0(:), y0(:), z0(:)];
start_pos = X_init(end,:)';
start_att = [0;0;0]*pi/180;
V0 = 10;

len0 = sqrt(sum((X_init(2:end,:)-X_init(1:end-1,:)).^2, 2) );
t_init = [reverse(-cumsum(reverse(len0./V0))); 0];

% energy_target = thermals(1,4:6)'+[0;0;-50]; %[400; 20; -150];
% energy_target = [info_target(1:2); -210];

%% SIMULATION VARIABLES
t0 = t_init(end);	% Initial time
tf = 2400;			% Final time
lookahead = 1;		% Lookahead horizon
t_plan = 5*lookahead;	% How long to plan ahead based on current wind estimate
replan = 3*lookahead;	% How frequently should the plan be recalculated
dt = 0.05;			% Time step for planning
ntf = 3;			% Number of path estimates

sim_dt = 0.01;	% Simulation timestep

n_xgrid = 8; n_ygrid = 8; n_zgrid = 5;
current_t = t0;

%% Turbulence - only generate if not already
% if ~exist('full_turb', 'var')
% 	fprintf(1, '\nGenerating turbulence...\n');
% 	[turb_uvw, turb_pqr] = ...
% 		dryden_cont(-z0(end), V0, b, u_20, (tf-t0)/sim_dt, sim_dt);
% 	full_turb = [turb_uvw'; turb_pqr'];
% 	fprintf(1, 'Done.\n');
% end

%% AUTOMATIC VARIABLES
X0 = zeros(12,1); X0(7:12) = [start_att; start_pos];

% --- Pitch rate offset cubic --- %
V_normal = 13; V_max = 40;
matmat = [V_normal^3, V_normal^2, V_normal, 1, 0;
		3*V_max^2, 2*V_max,     1   , 0, 0;
		  V_max.^3  ,  V_max.^2 ,  V_max  , 1,-1;
		  V_stall.^3, V_stall.^2,  V_stall, 1, 1];
pitch_coeff = rref(matmat); pitch_coeff = pitch_coeff(:,end)';

% --- Target flight vector --- %
info_target = start_pos;

n_replan = floor((tf-t0)/replan);		% Total number of replans
replan_points = round(replan/dt);		% Number of (long timestep) points

E0 = m*g*-start_pos(3) + 0.5*m*V0*V0;
current_pos = start_pos;	current_att = start_att;	current_V = V0;
old_pos = start_pos;		old_att = start_att;		old_V = V0;
current_E = E0;

obs_counter = 1;
t_full = t0:dt:tf-dt;
att_full = zeros(3, (tf-t0)/dt);
pos_full = zeros(3, (tf-t0)/dt);
V_full = zeros(1, (tf-t0)/dt);
X_train_full = zeros(max_obs, 3, n_replan);
W_train_full = zeros(max_obs, 3, n_replan);
t_train_full = zeros(max_obs, n_replan);
uncertainty = zeros(1, n_replan);

W_init = W_actual(X_init', t_init')' + wind_noise*randn(size(X_init));

% Training sets are initially just the manual observations. The training
% sets are updated as new data is recorded.
X_train = X_init;
W_train = W_init;
t_train = t_init;

% Get first energy target
[max_lift, energy_index] = min(W_train(:,3));
max_lift_old = max_lift;
energy_target = [X_train(energy_index, [1,2])'; current_pos(3)-10];
target_pos = energy_target;

% Grid to display current wind estimate
% [x_grid, y_grid, z_grid] = meshgrid( ...
% 	linspace(min(X_train(:,1)), max(X_train(:,1)), n_xgrid), ...
% 	linspace(min(X_train(:,2)), max(X_train(:,2)), n_ygrid), ...
% 	linspace(min(X_train(:,3)), max(X_train(:,3)), n_zgrid));
[x_grid, y_grid, z_grid] = meshgrid( ...
	linspace(x_limits(1),  x_limits(2), n_xgrid), ...
	linspace(y_limits(1),  y_limits(2), n_ygrid), ...
	linspace(z_limits(1),  z_limits(2), n_zgrid));
box_lims = [x_limits; y_limits; z_limits];


% Get first inv(K) matrix

if GPt
	Ki = GPt_predict(X_train, t_train, W_train, [], [], cov_funs{1}, loghyper);
	[WW_grid, VW_grid] = GPt_predict(X_train, t_train, W_train, ...
		[x_grid(:), y_grid(:), z_grid(:)], t0*ones(numel(x_grid), 1), ...
		cov_funs{1}, loghyper, Ki);
else
	Ki = GP_predict(X_train, W_train, [], cov_funs{1}, loghyper);
	[WW_grid, VW_grid] = GP_predict(X_train, W_train, ...
		[x_grid(:), y_grid(:), z_grid(:)], cov_funs{1}, loghyper, Ki);
end


U_grid = reshape(WW_grid(:,1), size(x_grid));
V_grid = reshape(WW_grid(:,2), size(x_grid));
W_grid = reshape(WW_grid(:,3), size(x_grid));
VW_grid = reshape(VW_grid, size(x_grid));
scale = sqrt(max(sum(WW_grid.^2, 2)))/max_wind*0.4;

%% Plot original wind field
figure(5); clf;
W_true_grid = W_actual([x_grid(:), y_grid(:), z_grid(:)]', ...
	t0*ones(1, numel(x_grid)));
set(gca, 'Zdir', 'reverse'); set(gca, 'Ydir', 'reverse'); daspect([1,1,1]);
h_W_true = coneplot(x_grid, y_grid, z_grid, ...
	reshape(W_true_grid(1,:)', size(x_grid)), ...
	reshape(W_true_grid(2,:)', size(x_grid)), ...
	reshape(W_true_grid(3,:)', size(x_grid)), ...
	x_grid, y_grid, z_grid, 0.4);
view(3);
set(h_W_true, 'EdgeColor', 'none');

%% FIGURE SETUP
h_fig2 = figure(2); clf; opos = get(h_fig2, 'Position');
opos(1:2) = min(opos(1:2), [1920-fig_width, 1100-fig_height]);
set(h_fig2, 'Position', [opos(1:2), fig_width, fig_height]);
h_fig = figure(1); clf;  opos = get(h_fig, 'Position');
opos(1:2) = min(opos(1:2), [1920-fig_width, 1100-fig_height]);
set(h_fig, 'Position', [opos(1:2), fig_width, fig_height]);
set(gca, 'Zdir', 'reverse'); set(gca, 'Ydir', 'reverse'); axis equal;
% h_plane = SBXC_handle(X0, PLANE_AERO);
view(3); hold on;
colours = [1 0.6 0.6; 0.6 1 0.6; 0.6 0.6 1; 1 0.6 1; 1 1 0.6; 0.6 1 1];
xlabel('X'); ylabel('Y'); zlabel('Z');

h_box = plot3(x_limits([1,2,2,1,1,1,2,2,1,1,1,1,2,2,2,2]), ...
			  y_limits([1,1,1,1,1,2,2,2,2,2,2,1,1,2,2,1]), ...
			  z_limits([1,1,2,2,1,1,1,2,2,1,2,2,2,2,1,1]), ...
			  'Color', [.8, .8, .8], 'LineStyle', '--');
h_target = plot3(target_pos(1), target_pos(2), target_pos(3), 'ro', ...
	'Color', [0.8,0,0], 'MarkerSize', 10, 'LineWidth', 2);
h_start = plot3(start_pos(1), start_pos(2), start_pos(3), 'g^', 'Color',...
	 [0,.5,0],'MarkerSize', 10, 'LineWidth', 2);
h_energy_target = plot3(energy_target(1), energy_target(2), ...
	energy_target(3), 'x','Color', [0.5,0.5,1], 'MarkerSize', 10);
h_info_target = plot3(info_target(1), info_target(2), ...
	info_target(3), 'bx','Color', [1,0.5,0], 'MarkerSize', 10);

h_train = plot3(X_train(:,1), X_train(:,2), X_train(:,3), 'r+', 'MarkerSize', 3);

h_windcone = coneplot(x_grid, y_grid, z_grid, U_grid, V_grid, W_grid, ...
	x_grid, y_grid, z_grid, scale, VW_grid);

set(h_windcone, 'EdgeColor', 'none')
set(gca, 'XLim', x_limits+[-100 100], 'YLim', y_limits+[-100 100], 'ZLim', z_limits+[-20 20]);

h_waitbar = waitbar(0, 'Wait, Peasant.');

%%

for ii = 1:n_replan
	
	if ~mod(ii, 5) && retrain && (ii>10)
		opt = optimset('Display', 'iter', 'GradObj', 'on', 'TolX', 1e-7,...
			'MaxIter', 100);
		if GPt
			[loghyper, fmin] = fminunc(@(hyper) ...
				GPt_likelihoodn(X_train, t_train, W_train(:,1), W_train(:,2), W_train(:,3), cov_funs, [hyper; loghyper0; sigma_hyper]), loghyper, opt);
			% fprintf(1, '\nRetrain: l = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', exp(loghyper));
			fprintf(1, ['\nRetrain: l = %0.5g, lt = %0.5g, ld = %0.5g sf = %0.5g, \n', ...
				'K = %0.5g, Wx = %0.5g, Wy = %0.5g, Wz = %0.5g, sn = %0.5g\n'], ...
				[exp(loghyper(1:4)), loghyper(5:8), exp(loghyper(9))]);
			
		else
			[loghyper, fmin] = fminunc(@(hyper) ...
				GP_likelihoodn(X_train, W_train(:,1), W_train(:,2), W_train(:,3), cov_funs, [hyper; loghyper0; sigma_hyper]), loghyper, opt);
			fprintf(1, '\nRetrain: l = %0.5g, sf = %0.5g, sn = %0.5g\n', exp(loghyper));
		end
		save_hypers(hyper_count,:) = [ii, loghyper];
		hyper_count = hyper_count+1;
		K_explore = 0.2/sqrt(exp(loghyper(end-1)));
	end
	
	if GPt
		Ki = GPt_predict(X_train, t_train, W_train, [], [], ...
			cov_funs{1}, loghyper);
		[WW_grid, VW_grid] = GPt_predict(X_train, t_train, W_train, ...
			[x_grid(:), y_grid(:), z_grid(:)], ...
			current_t*ones(numel(x_grid), 1), cov_funs{1}, loghyper, Ki);
	else
		Ki = GP_predict(X_train, W_train, [], cov_funs{1}, loghyper);
		[WW_grid, VW_grid] = GP_predict(X_train, W_train, ...
			[x_grid(:), y_grid(:), z_grid(:)], cov_funs{1}, loghyper, Ki);
	end
	
	U_grid = reshape(WW_grid(:,1), size(x_grid));
	V_grid = reshape(WW_grid(:,2), size(x_grid));
	W_grid = reshape(WW_grid(:,3), size(x_grid));
	VW_grid = reshape(VW_grid, size(x_grid));
	uncertainty(ii) = sum(VW_grid(:));
% 	delete(h_windcone); scale = sqrt(max(sum(WW_grid.^2, 2)))/max_wind*0.4;
% 	h_windcone = coneplot(x_grid, y_grid, z_grid, U_grid, V_grid, W_grid, ...
% 		x_grid, y_grid, z_grid, scale, VW_grid);
% 	set(h_windcone, 'EdgeColor', 'none')
	
	% Determine a target location based on energy and information estimate
	% from GP prediction
	[maxV, target_index] = max(VW_grid(:));
	info_target = [x_grid(target_index); y_grid(target_index); z_grid(target_index)];

	target_dist = sqrt(square_dist([x_grid(:), y_grid(:), z_grid(:)], current_pos'));
	[max_lift, energy_index] = min(W_grid(:)+VW_grid(:)+target_dist./2000);
	energy_target = ...
		[x_grid(energy_index); y_grid(energy_index); current_pos(3)+max_lift*50];
	%  Scale vertical height to soar in thermals but still have same height
	%  shear layer

	% Energy altitude required to get to info point
	Einfo = (sqrt(sum((info_target(1:2)-current_pos(1:2)).^2))/GR_approx + current_pos(3) - info_target(3))*m*g ;
	Einfo = Einfo*(max_lift < -0.5);
	
	target_pos = (Einfo > 0)*energy_target + (Einfo <= 0)*info_target;
	
	% Update target positions in plot
% 	set(h_target, 'XData', target_pos(1), 'YData', target_pos(2), ...
% 		'ZData', target_pos(3));
% 	set(h_energy_target, 'XData', energy_target(1), 'YData', energy_target(2), ...
% 		'ZData', energy_target(3));
% 	set(h_info_target, 'XData', info_target(1), 'YData', info_target(2), ...
% 		'ZData', info_target(3));
	
	% Calculate path based on target location and estimated wind field
	if GPt
		W_estimate = @(X_t, t_t) GPt_predict(X_train, t_train, W_train, ...
			X_t, t_t, cov_funs{1}, loghyper, Ki);
		J_estimate = @(X_t, t_t) dGPt_predict(X_train, t_train, W_train,...
			X_t, t_t, [1,2,3], cov_funs{3}, loghyper, Ki);
		
		[full_controls] = calculate_path_GP(current_pos, current_att, ...
			current_V, target_pos, W_estimate, J_estimate, t_plan, lookahead, dt, ntf, current_t);
	
	else
		W_estimate = @(X_test, t_t) GP_predict(X_train, W_train, X_test, cov_funs{1}, loghyper, Ki);
		J_estimate = @(X_test, t_t) dGP_predict(X_train, W_train, X_test, [1,2,3], cov_funs{3}, loghyper, Ki);
		
		[full_controls] = calculate_path_GP(current_pos, current_att, ...
			current_V, target_pos, W_estimate, J_estimate, t_plan, lookahead, dt, ntf);
	end
% 	axis tight;
	
% 	if movie_on; M((ii-1)*2 + 1) = getframe(h_fig, [0, 0, 640, 480]); end;
	
	% Execute the commands generated by the planner
	
	% With turbulence!
	
	[pos_path, att_path, V_path] = wind_sim3_control(current_pos, ...
		current_att, current_V, sim_dt, ...
		full_controls(:,1:round(replan/lookahead)), lookahead, ...
		W_actual, ntf);
	
	% Plot path and aircraft
% 	plot3(pos_path(1,:), pos_path(2,:), pos_path(3,:), 'k-');
% 	Cr = eye(4); Cr(1:3, 1:3) = calc_Ceb(att_path(:, end));	% Rotation	
% 	Ct = eye(4); Ct(1:3,4) = pos_path(:, end);			% Translation
% 	set(h_plane, 'Matrix', Ct*Cr)
	
	% New current position
	current_pos = pos_path(:,end);
	current_att = att_path(:,end);
	current_V = V_path(end);
	current_E = m*g*-current_pos(3) + 0.5*m*current_V.^2;
		
	% Fill in overall path
	pos_full(:,(ii-1)*replan_points+1:ii*replan_points) = pos_path(:,1:dt/sim_dt:end);
	att_full(:,(ii-1)*replan_points+1:ii*replan_points) = att_path(:,1:dt/sim_dt:end);
	  V_full(:,(ii-1)*replan_points+1:ii*replan_points) = V_path(:,1:dt/sim_dt:end);
	
	% Create new observation set
	n_obs = size(X_train, 1);
	obs_index = round((0:1/(freq_obs*replan):1)*size(pos_path, 2));
	obs_index = obs_index(2:end);
	new_X = pos_path(:,obs_index)';
	new_t = current_t + obs_index'*sim_dt;
	new_W = W_actual(new_X', new_t')' + wind_noise*randn(size(new_X));
	
	X_train_full(1:n_obs, :, ii) = X_train;
	W_train_full(1:n_obs, :, ii) = W_train;
	t_train_full(1:n_obs, ii) = t_train;
		
	% If the observation set is not full yet, add the new data
	if n_obs < max_obs
		if n_obs + size(new_X, 1) > max_obs
			X_train = [X_train; new_X(1:max_obs-n_obs,:)];
			W_train = [W_train; new_W(1:max_obs-n_obs,:)];
			t_train = [t_train; new_t(1:max_obs-n_obs)];
		else
			X_train = [X_train; new_X];
			W_train = [W_train; new_W];
			t_train = [t_train; new_t];
		end
		
% 		max_obs = max(size(X_train, 1), max_obs);
		
	% Otherwise, calculate the new training set
	else
		% Check if any points are outside useful region
		out_box = (X_train < ones(n_obs,1)*(box_lims(:,1)'-2*exp(loghyper(1))))|...
			(X_train > ones(n_obs,1)*(box_lims(:,2)'+2*exp(loghyper(1))));
		
		% Get distances between points (upper triangular due to symmetry)
		obs_d = triu(square_dist(X_train, X_train), 1);
		min_d = min(obs_d + tril(ones(n_obs))*1e6, [], 2);
		
		min_d = min_d.*(~any(out_box,2));
		
		% INCLUDE TIME OFFSET (divide by time since now)
		if GPt
			min_d = min_d./(current_t - t_train);
		end
		[sort_d, sort_i] = sort(min_d);
		for j = 1:size(new_X, 1)
			X_train(sort_i(j),:) = new_X(j,:);
			W_train(sort_i(j),:) = new_W(j,:);
			t_train(sort_i(j),:) = new_t(j,:);
		end
	end
% 	set(h_train, 'Xdata', X_train(:,1), 'YData', X_train(:,2), ...
% 		'ZData', X_train(:,3));
	
	
% 	if movie_on; M((ii-1)*2 + 2) = getframe(h_fig, [0, 0, 640, 480]); end;
	if movie_on
		M(ii) = getframe(h_fig, [0, 0, 640, 480]);
		set(0, 'CurrentFigure', h_fig2);
		M2(ii) = plot_variance(x_grid, y_grid, z_grid, VW_grid, [X_init', ...
			pos_full(:, 1:ii*replan_points)], box_lims, info_target, 5);
		set(0, 'CurrentFigure', h_fig);
	end
	
	current_t = current_t+replan;
	waitbar(ii/n_replan, h_waitbar);
end

%%
pos_full = pos_full(:,1:ii*replan_points);
att_full = att_full(:,1:ii*replan_points);
V_full	 = V_full(1:ii*replan_points);
t_full	 = t_full(1:ii*replan_points);

save_hypers = save_hypers(1:hyper_count-1,:);

%%
M(numel(M)+1) = getframe(h_fig, [0, 0, fig_width, fig_height]);
movie2avi(M, 'movies\soaring_OVERWRITE.avi', 'fps', 10, ...
'Compression','none')
movie2avi(M2, 'movies\soaring_OVERWRITE_isosurf.avi', 'fps', 10, ...
'Compression','none')

%%
figure(3); clf; hold on;
Ek = 0.5*m*V_full.^2;
Eks = Ek - Ek(1);
Ep = m*g*-pos_full(3,:);
Eps = Ep - Ep(1);
Et = Eks + Eps;
plot(t_full, Eks, 'b--', t_full, Eps, 'r-.');
hold on; plot(t_full, Et, 'Color', [0, .6, 0]);
grid on;
xlabel('Time, \it{t}\rm (s)');
ylabel('Energy change, \it\DeltaE \rm(J)');
h_l = legend('Kinetic energy', 'Potential energy', 'Total energy');
set(h_l, 'location', 'best');

%%
figure(4); clf;
plot(t0:replan:(n_replan-1)*replan, uncertainty./numel(VW_grid), 'b-');
grid on;
xlabel('Time, \it{t}\rm (s)');
ylabel('Mean field variance,  \it{mean}\rm(\sigma_i^2), (m^2/s^2)');
temp = get(gca, 'Ylim'); set(gca, 'YLim', [0, temp(2)]);


%% overplot
figure(5); clf;
[axx, ha1, ha2] = plotyy(t_full, Eps, t0:replan:(n_replan-1)*replan, uncertainty./numel(VW_grid));
set(ha1, 'linestyle', '-.', 'color', 'r');
set(axx(1), 'nextplot', 'add', 'ycolor', 'k');
plot(axx(1), t_full, Eks, 'b--'); plot(axx(1), t_full, Et, 'Color', [0, .6, 0]);

set(ha2, 'color', [0.9, 0.5, 0.1], 'linewidth', 2.0); set(axx(2), 'ycolor', [0.9, 0.5, 0.1]);

set(get(axx(1),'Xlabel'),'String','Time, \it{t}\rm (s)');
set(get(axx(1),'Ylabel'),'String','Energy change, \it\DeltaE \rm(J)') 
set(get(axx(2),'Ylabel'),'String','Mean field variance,  \it{mean}\rm(\sigma_i^2), (m^2/s^2)')

% axis tight
grid on

%% SAVE DATA
save Data\TEMP.mat
fprintf(1, '\nData saved to Data\\TEMP.mat \n');


%% SCRAP
% Old observationn set calculator
% 		for p = 1:size(new_X, 2)
% 			X_train(obs_counter,:) = new_X(p,:);
% 			W_train(obs_counter,:) = new_W(p,:);
% 			obs_counter = obs_counter+1;
% 			obs_counter = rem(obs_counter, max_obs)+1;
% 		end
%
% Sum distance calculator (no good)
% 		obs_d = square_dist(X_train, X_train);
% 		sum_d = sum(obs_d, 1);
% 		[sort_d, sort_i] = sort(sum_d);
% 		for j = 1:size(new_X, 1)
% 			X_train(sort_i(j),:) = new_X(j,:);
% 		end
