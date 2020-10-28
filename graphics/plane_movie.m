function [M] = plane_movie(plane_aero, frame_skip, X, U, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		plane_movie
%
% PURPOSE:		plot a movie of an aircraft flight
%               
% SYNTAX:		[M] = plane_movie(plane_aero, X, U, frame_skip)
%
% INPUTS:		plane_aero	- Aircraft definition structure
%				X			- Matrix of state vectors
%				U			- Corresponding control vectors
%				frame_skip	- Number of frames to skip
%
% OUTPUTS:		M			- MATLAB Movie file
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2005
%
% MODIFIED:     October 2007
%
% See also:		aero_plot, strip_plot
%--------------------------------------------------------------------------

set(gcf, 'Name', 'Attitude video');
set(gcf, 'position', [100 100 640 480]);

lag_frame = 5;			% Set the number of frames that the camera lags
lag_vec = [-12; 7; -1];	% Set the lag vector of the camera
offset = [0; 0; 4];	% Offset target to allow display

% Check number of planes to animate and whether wind field has been
% specified
if nargin < 4
	error('wing_model:plane_movie:nargin', 'Not enough input arguments');
elseif nargin == 4
	col = 'b-';
	nframes = size(X, 2);
else
	% nplane is the number of planes present, and iswind flags whether a
	% wind field has been specified.
	col = varargin{1};
	nplane = (nargin - 2)/3;
	if nplane ~= round(nplane)
		error('wing_model:plane_movie:typeargin',...
			['Input must be in the form: plane_movie(aero, fskip, X0, ',...
			'U0, col0, X1, U1, col1, ..., XYZ, UVW, scale']);
	end
	
	% Determine if wind is specified by checking last argument (if it is a
	% string then it is the colour to draw the plane, otherwise it must be
	% a wind field scale).
	iswind = ~ischar(varargin{end});
	if iswind
		XYZ = varargin{end-2};
		UVW = varargin{end-1};
		scale = varargin{end};
	end
	
	nplane = nplane - iswind;
	nframes = zeros(1,nplane);
	nframes(1) = size(X, 2);
	
	% If there are multiple planes, sort out the arguments
	if nplane > 1
		X_extra   = cell(1, nplane-1);
		U_extra   = cell(1, nplane-1);
		col_extra = cell(1, nplane-1);
		for i = 1:nplane-1
			X_extra{i}  = varargin{2+(i-1)*3};
			U_extra{i}  = varargin{3+(i-1)*3};
			col_extra{i}= varargin{4+(i-1)*3};
			nframes(i+1) = size(X_extra{i}, 2);
		end
	end

end

h_axis = zeros(1,3);

% ------------------------------------------------------------------ %
% Initialise axes
%	h_axis(1) is flight view
%	h_axis(2) is trajectory,
%	h_axis(3) is artificial horizon
%	h_axis(4) is heading indicator
%	h_axis(5) is slip/skid indicator
%		- All axes are children of the current figure
% ------------------------------------------------------------------ %

% ------------------------------------------------------------------ %
% Intialise image in main axes (h_axis(1))
% ------------------------------------------------------------------ %
h_axis(1) = axes('Position', [0, 0, 1, 1]);
hold on; grid on;
plot3(X(10,:), X(11,:), X(12,:), col);
axis image; view(3); axis(axis + [-10, 10, -10, 10, -10, 10])
% set(gca,'YDir','reverse'); 	set(gca,'ZDir','reverse');
h_plane = zeros(1, nplane);
h_plane(1) = SBXC_handle(zeros(12,1), plane_aero);
for i = 1:nplane-1
	h_plane(i+1) = SBXC_handle(X_extra{i}(:,1), plane_aero);
end

if any(X(1:3,1) ~= 0)
	[FF, MM, h_arrows] = strip_method(X(:,1), U(:,1), 1);
end
if iswind
	h_quiver2 = quiver3(XYZ{1}, XYZ{2}, XYZ{3}, UVW{1}, UVW{2}, UVW{3}, scale);
	set(h_quiver2, 'color', [0.6 0 0]);
end


set(h_axis(1), 'CameraPosition', X(10:12,1)+lag_vec)
set(h_axis(1), 'CameraTarget', X(10:12,1)+offset)
set(h_axis(1), 'CameraUpVector', [0;0;-1])
set(h_axis(1), 'CameraViewAngle', 60)


% ------------------------------------------------------------------ %
% Intialise trajectory view (h_axis(2))
% ------------------------------------------------------------------ %
h_axis(2) = axes('Position', [0, 0, 0.4, 0.4]);
hold on; grid on;
plot3(X(10,:), X(11,:), X(12,:), col);
axis image; view(3); axis(axis + [-10, 10, -10, 10, -10, 10])
set(gca,'YDir','reverse'); set(gca,'ZDir','reverse');
h_dot = zeros(1, nplane);
h_dot(1) = plot3(X(10,1), X(11,1), X(12,1), 'r.', 'MarkerSize', 15);
for i = 1:nplane-1
	plot3(X_extra{i}(10,:), X_extra{i}(11,:), X_extra{i}(12,:), col_extra{i});
	h_dot(i+1) = plot3(X_extra{i}(10,1), X_extra{i}(11,1), ...
		X_extra{i}(12,1), 'r.', 'MarkerSize', 15);
end
if iswind
	h_quiver2 = quiver3(XYZ{1}, XYZ{2}, XYZ{3}, UVW{1}, UVW{2}, UVW{3}, scale*5);
	set(h_quiver2, 'color', [0.6 0 0]);
end

% ------------------------------------------------------------------ %
% Initialise artificial horizon
% ------------------------------------------------------------------ %
horiz_size = 0.25;
h_axis(3) = axes('Position', [1-horiz_size, 0.05, horiz_size, horiz_size]);
axis([-1 1 -1 1]); axis equal; axis off; hold on;
horizon_handle = artificial_horizon(X(7,1), X(8,1));

% ------------------------------------------------------------------ %
% Initialise heading indicator
% ------------------------------------------------------------------ %
h_axis(4) = axes('Position', [1-horiz_size, 0, horiz_size, 0.05]);
axis([-1 1 -1 1]); hold on;
set(h_axis(4), 'Color', 'k', 'XTick', [], 'YTick', []);
heading_handle = heading_indicator(X(9,1));

% ------------------------------------------------------------------ %
% Initialise slip/skid ball
% ------------------------------------------------------------------ %
h_axis(5) = axes('Position', [0.4, 0.01, 0.2, 0.05]);
beta_max = 10;
axis([-beta_max beta_max -0.5 0.5]); hold on;
set(h_axis(5), 'XTick', [], 'YTick', []);
plot([-1 1 1 -1]*beta_max, [-1 -1 1 1]*0.5, '-k', 'LineWidth', 2);
plot([-0.5 -0.5], [-0.5 0.5], 'k:', [0.5 0.5], [-0.5 0.5], 'k:', ...
	[-5 -5], [-0.5 0.5], 'k:', [5 5], [-0.5 0.5], 'k:');
slip_handle = slip_indicator(X(:,1), beta_max);

set(gcf, 'CurrentAxes', h_axis(1));
M(1) = getframe(gcf);
k=2;
tstart = tic;
for ii = 2:(frame_skip+1):max(nframes)
	
	jframe = ii.*(ii < nframes) + nframes.*(ii >= nframes);
	i = jframe(1);
	Cr = eye(4);
	Cr(1:3, 1:3) = calc_Ceb(X(7:9,i));		% Rotation transformation	
	Ct = eye(4);
	Ct(1:3,4) = X(10:12, i);		% Translate to actual position
	
	% Slip indicator
	set(gcf, 'CurrentAxes', h_axis(5));
	delete(slip_handle);
	slip_handle = slip_indicator(X(:,i), beta_max);
	toc(tstart);
	
	% Heading indicator
	set(gcf, 'CurrentAxes', h_axis(4));
	delete(heading_handle);
	heading_handle = heading_indicator(X(9,i));
	toc(tstart);
	
	% Artificial horizon
	set(gcf, 'CurrentAxes', h_axis(3));
	delete(horizon_handle);
	horizon_handle = artificial_horizon(X(7,i), X(8,i));
	toc(tstart);
	
	% Trajectory view
	set(gcf, 'CurrentAxes', h_axis(2));	
	delete(h_dot);
	h_dot(1) = plot3(X(10,i), X(11,i), X(12,i), 'r.', 'MarkerSize', 15);
	for j = 1:nplane-1
		jf = jframe(j+1);
		h_dot(j+1) = plot3(X_extra{j}(10,jf), X_extra{j}(11,jf), ...
			X_extra{j}(12,jf), 'r.', 'MarkerSize', 15);
	end
	toc(tstart);
	
	% Main flight view
	set(gcf, 'CurrentAxes', h_axis(1));
	if exist('h_arrows', 'var')
		if h_arrows{1} ~= 0; delete(h_arrows{1}); end
		if h_arrows{2} ~= 0; delete(h_arrows{2}); end
    end
	
%     CC = makehgtform('xrotate', X(7,i), 'yrotate', X(8,i), ...
% 		'zrotate', X(9,i), 'translate', X(10,i), X(11,i), X(12,i));
	set(h_plane(1), 'Matrix', Ct*Cr);
% 	toc(tstart);
	

	if i > lag_frame
% 		Ceb = calc_Ceb(X(7:9,i-lag_frame));
		cam_pos = X(10:12,i) +lag_vec;% + calc_Ceb(X(7:9,i))*lag_vec;
		set(h_axis(1), 'CameraPosition', cam_pos);
		set(h_axis(1), 'CameraTarget', X(10:12,i)+offset)
% 		set(h_axis(1), 'CameraUpVector', Ceb*[0;0;-1])
	else
		set(h_axis(1), 'CameraTarget', X(10:12,i)+offset)
		set(h_axis(1), 'CameraPosition', X(10:12,i)+lag_vec);
	end
	
	if any(X(1:3,i) ~= 0)
		[FF, MM, h_arrows] = strip_method(X(:,i), U(:,i), 1);
	end
	
	for j = 1:nplane-1
		jf = jframe(j+1);
		Cr = eye(4);
		Cr(1:3, 1:3) = calc_Ceb(X_extra{j}(7:9,jf));		% Rotation transformation	
		Ct = eye(4);
		Ct(1:3,4) = X_extra{j}(10:12, jf);		% Translate to actual position
		set(h_plane(j+1), 'Matrix', Ct*Cr);
	end
	
	drawnow;
	M(k) = getframe(gcf);
% 	toc(tstart);
	k = k + 1;
end   

% M = M(1:k-2);
% h_ax_fake = axes('Position', [0 0 1 1], 'Visible', 'off');
set(gcf, 'CurrentAxes', h_axis(1));
% close(mov);
% movie(M, [1:length(M)], 10)


% NOTES: TO MAKE AVI MOVIES
% movie2avi(M, 'movies\MOVIE_NAME.avi', 'fps', 10, 'Compression','none')
% To release all movies for deletions and to prevent sharing violations,
% clear mex