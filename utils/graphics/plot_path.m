function plot_path(X, plane_aero, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		plot_path
%
% PURPOSE:		Plot aircraft path with cycles transitions if specified
%
% SYNTAX:		plot_path(X, aero)
%				plot_path(X, aero, tx, mode)
%				plot_path(X, aero, tx, mode, ls, col)
%				plot_path(X, aero, tx, mode, plane_int, ls, col)
%
% INPUTS:		X		- [n×m] State vector - n states, m points in time
%				aero	- aerodynamic structure of plane (see SBXC_def)
%				tx		- [1×m] Timing vector for X
%				mode	- dynamic soaring mode change times
%				plane_int time interval to plot plane pose samples
%				ls		- line style (standard styles, '-', '-.', etc.)
%				col		- colour (standard colors 'b', 'g' etc or [1x3] RGB
%
% OUTPUTS:		(graphical - current axis)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2008
%
% MODIFIED:     July 2008
%
% See also:		plot_states, SBXC_def
%--------------------------------------------------------------------------
modes_flag = 0;
ls = '-'; col = 'b';
fsize = 8;

if nargin > 3; modes_flag = 1; end;

if nargin == 4
	tx = varargin{1};
	modetimes = varargin{2};
	plane_int = tx(end)+1;
elseif nargin == 5
	tx = varargin{1};
	modetimes = varargin{2};
	plane_int = varargin{3};
elseif nargin == 6
	tx = varargin{1};
	modetimes = varargin{2};
	plane_int = tx(end)+1;
	ls = varargin{3};
	col = varargin{4};
elseif nargin == 7
	tx = varargin{1};
	modetimes = varargin{2};
	plane_int = varargin{3};
	ls = varargin{4};
	col = varargin{5};
end

SBXC_handle(X(:,end), plane_aero, 1);
set(gca, 'Xdir', 'normal', 'YDir', 'normal', 'ZDir', 'normal')
set(gca, 'FontSize', fsize)

plot3(X(10,:), X(11,:), X(12,:), ls, 'Color', col); hold on;
plot3(X(10,1), X(11,1), X(12,1), 'g^');
plot3(X(10,end), X(11,end), X(12,end), 'ro');

if modes_flag
	
	% Plot points for cycle changes
	ncycles = size(modetimes, 1);
	cycle_times = zeros(1, ncycles);
	for i = 1:ncycles
		cycle_times(i) = find(abs(tx - modetimes(i,1)) < 1e-5);
	end
	plot3(X(10, cycle_times), X(11, cycle_times), X(12, cycle_times), ...
		'.', 'Color', [1, .5, 0])
	
	% Plot points for mode changes
	modechange = modetimes(:,2:end);
	nchanges = numel(modechange);
	change_times = zeros(1, nchanges);
	for i = 1:nchanges
		change_times(i) = find(abs(tx - modechange(i)) < 1e-5);
	end
	plot3(X(10, change_times), X(11, change_times), X(12, change_times), ...
		'.', 'Color', [0, .5, .5])
	
	n_intervals = floor(tx(end)/plane_int);
	for i = 0:n_intervals
		t_current = find(abs(tx-i*plane_int) < 1e-8, 1);
		SBXC_handle(X(:,t_current), plane_aero, 1);
	end
end
set(gca, 'Xdir', 'normal', 'YDir', 'reverse', 'ZDir', 'reverse')
axis equal;
axis([min(X(10,:)), max(X(10,:)), min(X(11,:)), max(X(11,:)), ...
	min(X(12,:)), max(X(12,:))] + 10*[-1, 1, -1, 1, -1, 1]);