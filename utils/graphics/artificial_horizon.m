function horizon_handle = artificial_horizon(bank, pitch)

pitch_max = 50*pi/180;

pitch_ratio = pitch/pitch_max;

sky_blue = [0, 0.5, 1];
ground_brown = [180, 90, 10]./255;

% if abs(pitch_ratio) >= 1
% 	theta_start = bank - pi*(sign(pitch)+1)/2 + 5*pi/180;
% 	theta_stop  = bank - pi*(sign(pitch)+1)/2 - 5*pi/180;
% else
	theta_start = bank - real(asin(pitch_ratio));
	theta_stop  = pi + bank + real(asin(pitch_ratio));
% end

	
% Plot horizon
hor_handle(1) = plot_semi(theta_start, theta_stop, 1, sky_blue);
hor_handle(2) = plot_semi(theta_stop, theta_start, 1, ground_brown);
hor_handle(3) = plot([cos(theta_stop), cos(theta_start)], ...
	[sin(theta_stop), sin(theta_start)], 'Color', [0 0 0], 'LineWidth', 2);

% Plot wings
wing_handle(1) = plot([-.7 -.1 0 .1 .7], [0 0 -.1 0 0], 'y', 'LineWidth', 2);
wing_handle(2) = plot([0], [0], 'y+');

% Plot pitch indicators
min_marker = -30*pi/180;
max_marker = 30*pi/180;
int_marker = 10*pi/180;
len = 0.4;
k = 1;
pitchl_handle = zeros(1, 2*(max_marker-min_marker)/int_marker);

% Minor pitch lines (5 degree increments)
for pitch_line = min_marker:int_marker:max_marker
	d = pitch_line/pitch_max;
	pitchl_handle(k) = plot(...
		[-d*sin(bank) - len/2*cos(bank), -d*sin(bank) + len/2*cos(bank)], ...
		[d*cos(bank) - len/2*sin(bank), d*cos(bank) + len/2*sin(bank)], ...
		'Color', [1 1 1], 'LineWidth', 1);
	
	k = k+1;
	d = d - int_marker/2/pitch_max;
	pitchl_handle(k) = plot(...
		[-d*sin(bank) - len/4*cos(bank), -d*sin(bank) + len/4*cos(bank)], ...
		[d*cos(bank) - len/4*sin(bank), d*cos(bank) + len/4*sin(bank)], ...
		'Color', [1 1 1], 'LineWidth', 1);
	
	k = k+1;
end
	
horizon_handle = [hor_handle, wing_handle, pitchl_handle];

% 	horizon_handle = plot([cos(theta_stop), cos(theta_start)], ...
% 	[sin(theta_stop), sin(theta_start)]);
