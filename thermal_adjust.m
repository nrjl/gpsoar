function t_param = thermal_adjust(t_param, V_max, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		thermal_adjust
%
% PURPOSE:		Adjust thermal group to meet total lift volume requirement.
%				Can specify either total lifting colume requirement V_lift
%				or the parameters for an equivalent single thermal. If
%				current lifting volume is too much, remove weakest thermal
%				and repeat, otherwise increase speed of last hermal to
%				match.
%
% SYNTAX:		t_param = thermal_adjust(t_param, V_lift)
%				t_param = thermal_adjust(t_param, V_max, R_max)
%
% INPUTS:		t_param	- Thermal parameters
%				V_lift	- Total lifting volume required
%				V_max	- Equivalent max vertical speed
%				V_max	- Equivalent radius
%
% OUTPUTS:		t_param	- Adjusted thermal parameters
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2010
%
% MODIFIED:     November 2010
%
% See also:		
%--------------------------------------------------------------------------
if nargin == 2
	target_lift_volume = V_max;
else
	target_lift_volume = 4*varargin{1}^2*V_max/pi;
end
n_thermals = size(t_param, 1);

if n_thermals == 1
	t_param(1) = target_lift_volume*pi/4/(t_param(2).^2);
	return
end

lift_volume = 4/pi*t_param(:,2).^2.*t_param(:,1);
actual_lift_volume = sum(lift_volume);


if actual_lift_volume > target_lift_volume
	[maxlift, bigtherm] = max(lift_volume);
	if (actual_lift_volume - maxlift) > target_lift_volume
		t_param(bigtherm,:) = [];
		t_param = thermal_adjust(t_param, target_lift_volume);
	else
		% Calculate the actual lift required by the (old) largest thermal
		required_lift = target_lift_volume - (actual_lift_volume-maxlift);
		t_param(bigtherm, 1) = required_lift*pi/4/(t_param(bigtherm,2).^2);
	end
else
	required_lift = target_lift_volume - ...
		sum(4/pi*t_param(1:(n_thermals-1),2).^2.*t_param(1:(n_thermals-1),1));
	t_param(n_thermals, 1) = required_lift*pi/4/(t_param(n_thermals,2).^2);
end
