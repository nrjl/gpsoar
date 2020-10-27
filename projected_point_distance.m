function d = projected_point_distance(x, t, y, x_star, t_star)
%--------------------------------------------------------------------------
%
% FUNCTION:		projected_point_distance
%
% PURPOSE:		Find the distance from a projected vector to a target
%				point
%               
% SYNTAX:		d = projected_point_distance(x, t, y, x_star, t_star)
%
% INPUTS:		x	- observation location(s) - [n×d]
%				t	- observation time(s)	- [n×1]
%				y	- vector observation(s) - [n×d]
%				x_star	- target point(s)	- [k×d]
%				t_star	- target time(s)	- [1×1]
%
% OUTPUTS:		d	- square distances [n×k]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		March 2010
%
% MODIFIED:     March 2009
%
% See also:		
%--------------------------------------------------------------------------

% Now can only predict at one time

% Calculate projected position (position of observation at target time)
% New_pos	= x + y*(t_star - t)
%			= (x - y*t) + y*t_star

projected_pos = permute(repmat(x + (t_star - diag(t))*y, [1, 1, size(x_star, 1)]), [1, 3, 2]);

target_pos = permute(repmat(x_star, [1, 1, size(x, 1)]), [3, 1, 2]);

d = sum((target_pos - projected_pos).^2, 3);


% OLD VERSION (NON - SYMMETRIC)
% start_pos = permute(repmat(x - diag(t)*y, [1, 1, size(x_star, 1)]), [1, 3, 2]);
% shift_pos = zeros(size(start_pos));
% 
% % Horrible hack to fill dimensions of shifted positions
% for ii = 1:size(x, 2)
% 	shift_pos(:,:,ii) = y(:,ii)*t_star';
% end
% 
% projected_pos = start_pos + shift_pos;
% target_pos = permute(repmat(x_star, [1, 1, size(x, 1)]), [3, 1, 2]);
% 
% d = sum((target_pos - projected_pos).^2, 3);