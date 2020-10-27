function V_wg_e = wind_field(locs, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		wind_field
%
% PURPOSE:		Calculate a wind field at a set of locations in space
%               
% INPUTS:		locs	- [3 x n] matrix of earth frame locations
%
% OUTPUTS:		V_wg_e	- Velocity of wind with respect to ground in earth
%						  frame coordinates (NED)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		September 2007
%
% MODIFIED:		September 2007
%
% See also:		blade_element, strip_forces
%--------------------------------------------------------------------------
global psi_w			% Wind bearing angle

if nargin == 2
	V_wg_e = zeros(3, size(locs, 2));
	V_wg_e(3,:) = varargin{1}*ones(1, size(locs, 2));
else
	
	V_wg_e = zeros(size(locs));
	
 	V = linear_profile(-locs(3,:), 100, 115, 0, 15);
% 	V = cos_profile(-locs(3,:), 100, 120, 1, 15);
%	V = pohlhausen(-locs(3,:), 100, 12, 15, 1);
% 	[Wx, Wy] = linear_wind(locs(1,:), locs(2,:), 0, 0, [0 0.1; 0.0 0]);
% 
% 	V_wg_e(1,:) = Wx*cos(psi_w);
% 	V_wg_e(2,:) = Wx*sin(psi_w);
% 	V_wg_e(3,:) = -Wy;
	
	Wx = V*cos(psi_w);
	Wy = V*sin(psi_w);

	V_wg_e(1,:) = Wx;
	V_wg_e(2,:) = Wy;
	%V_wg_e(3,:) = zeros(1, size(locs,2));
% 	
% 
% % 	V_wg_e(3,:) = (cos(2*pi*locs(1,:)./20)-1)*-0.1;
% % 	V_wg_e(3,:) = -0.5+acos(cos(2*pi*locs(1,:)./80))/pi;

    
%     V_wg_e(1,:) = -0.5*locs(1,:);
end

