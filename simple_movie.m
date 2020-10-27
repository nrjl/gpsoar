function [M] = simple_movie(plane_aero, frame_skip, pos, att)
%--------------------------------------------------------------------------
%
% FUNCTION:		simple_movie
%
% PURPOSE:		plot a simple movie of an aircraft flight
%               
% SYNTAX:		[M] = simple_movie(plane_aero, frame_skip, pos, att)
%
% INPUTS:		plane_aero	- Aircraft definition structure
%				frame_skip	- Number of frames to skip
%
% OUTPUTS:		M			- MATLAB Movie file
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2009
%
% MODIFIED:     October 2009
%
% See also:		plane_movie
%--------------------------------------------------------------------------

clf;
set(gcf, 'Name', 'Attitude video');
set(gcf, 'position', [100 100 640 480]);

lag_frame = 5;			% Set the number of frames that the camera lags
lag_vec = [-12; 7; -1];	% Set the lag vector of the camera
offset = [0; 0; 4];	% Offset target to allow display

% ------------------------------------------------------------------ %
% Intialise image in main axes (h_axis(1))
% ------------------------------------------------------------------ %
h_axis(1) = axes('Position', [0, 0, 1, 1]);
hold on; grid on;
plot3(pos(1,:), pos(2,:), pos(3,:), 'r-');
set(gca,'ZDir','reverse'); set(gca,'YDir','reverse'); 
axis image; view(3); axis(axis + [-10, 10, -10, 10, -10, 10])
h_plane = SBXC_handle(zeros(12,1), plane_aero);

set(h_axis(1), 'CameraPosition', pos(:,1)+lag_vec)
set(h_axis(1), 'CameraTarget', pos(:,1)+offset)
set(h_axis(1), 'CameraUpVector', [0;0;-1])
set(h_axis(1), 'CameraViewAngle', 60)

% ------------------------------------------------------------------ %
% Intialise trajectory view (h_axis(2))
% ------------------------------------------------------------------ %
h_axis(2) = axes('Position', [0, 0, 0.4, 0.4]);
hold on; grid on;
plot3(pos(1,:), pos(2,:), pos(3,:), 'r-');
set(gca,'ZDir','reverse'); set(gca,'YDir','reverse'); 
axis image; view(3); axis(axis + [-10, 10, -10, 10, -10, 10])
h_dot = plot3(pos(1,1), pos(2,1), pos(3,1), 'r.', 'MarkerSize', 15);

% ------------------------------------------------------------------ %
% Initialise artificial horizon
% ------------------------------------------------------------------ %
horiz_size = 0.25;
h_axis(3) = axes('Position', [1-horiz_size, 0.05, horiz_size, horiz_size]);
axis([-1 1 -1 1]); axis equal; axis off; hold on;
horizon_handle = artificial_horizon(att(1,1), att(2,1));

% ------------------------------------------------------------------ %
% Initialise heading indicator
% ------------------------------------------------------------------ %
h_axis(4) = axes('Position', [1-horiz_size, 0, horiz_size, 0.05]);
axis([-1 1 -1 1]); hold on;
set(h_axis(4), 'Color', 'k', 'XTick', [], 'YTick', []);
heading_handle = heading_indicator(att(3,1));

set(gcf, 'CurrentAxes', h_axis(1));
M(1) = getframe(gcf);
nframes = size(pos,2);
k=2;

for ii = 2:(frame_skip+1):nframes
	
	jframe = ii.*(ii < nframes) + nframes.*(ii >= nframes);
	i = jframe(1);
	Cr = eye(4);
	Cr(1:3, 1:3) = calc_Ceb(att(:,i));		% Rotation transformation	
	Ct = eye(4);
	Ct(1:3,4) = pos(:, i);		% Translate to actual position
	
	% Heading indicator
	set(gcf, 'CurrentAxes', h_axis(4));
	delete(heading_handle);
	heading_handle = heading_indicator(att(3,i));
	
	% Artificial horizon
	set(gcf, 'CurrentAxes', h_axis(3));
	delete(horizon_handle);
	horizon_handle = artificial_horizon(att(1,i), att(2,i));
	
	% Trajectory view
	set(gcf, 'CurrentAxes', h_axis(2));	
	delete(h_dot);
	h_dot(1) = plot3(pos(1,i), pos(2,i), pos(3,i), 'r.', 'MarkerSize', 15);
	
	% Main flight view
	set(gcf, 'CurrentAxes', h_axis(1));
	
% 	CC = makehgtform('xrotate', X(7,i), 'yrotate', X(8,i), ...
% 		'zrotate', X(9,i), 'translate', X(10,i), X(11,i), X(12,i));
	set(h_plane(1), 'Matrix', Ct*Cr);
	
	if i > lag_frame
		cam_pos = pos(:,i) +lag_vec;
		set(h_axis(1), 'CameraPosition', cam_pos);
		set(h_axis(1), 'CameraTarget', pos(:,i)+offset)
	else
		set(h_axis(1), 'CameraTarget', pos(:,i)+offset)
		set(h_axis(1), 'CameraPosition', pos(:,i)+lag_vec);
	end
	
	drawnow;
 	M(k) = getframe(gcf);
	k = k + 1;
end   

set(gcf, 'CurrentAxes', h_axis(1));

% NOTES: TO MAKE AVI MOVIES
% movie2avi(M, 'movies\MOVIE_NAME.avi', 'fps', 10, 'Compression','none')
% To release all movies for deletions and to prevent sharing violations,
% clear mex