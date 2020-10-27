function C = calc_Ceb(varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		calc_Ceb
%
% PURPOSE:		calculate the transformation matrix from body to Earth axes
%               
% SYNTAX:		C = calc_Ceb(phi, theta, psi)
%				C = calc_Ceb(angles)	where angles is a 1×3 vector
%
% INPUTS:		phi		- roll angle (rad)
%				theta	- pitch angle (rad)
%				psi		- heading angle (rad)
%				angles	- [phi, theta, psi]
%
% OUTPUTS:		axis_lim- axis limits
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		March 2007
%
% MODIFIED:     June 2007
%
% See also:		sbxc_full, plot_polars
%--------------------------------------------------------------------------

if nargin == 1
	angles = varargin{1};
	angles = angles(:);
	
	if max(size(angles) ~= [3, 1])
		error('blade_element:calc_Ceb:in_size', ...
			'Angle input must contain three angles')
	end
	
	phi = angles(1);
	theta = angles(2);
	psi = angles(3);
	
elseif nargin == 3
	phi = varargin{1};
	theta = varargin{2};
	psi = varargin{3};
else
	error('blade_element:calc_Ceb:n_inputs', ['Three angles must be ', ...
		'entered, either as a vector or separate values']);
end
	
	
C = [cos(theta)*cos(psi)	sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi)	cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
	cos(theta)*sin(psi)		sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi)	cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
	-sin(theta)				sin(phi)*cos(theta)								cos(phi)*cos(theta)	];
	
		
		



