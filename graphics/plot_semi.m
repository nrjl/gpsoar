function out_handle = plot_semi(theta_start, theta_stop, R, col)

if theta_stop < theta_start
	theta_stop = theta_stop+2*pi;
end

n_interval = round((theta_stop - theta_start)/pi*180)+1;
theta_vec = linspace(theta_start, theta_stop, n_interval);

x = R*cos(theta_vec);
y = R*sin(theta_vec);

out_handle = fill([x, x(1)], [y, y(1)], col);

