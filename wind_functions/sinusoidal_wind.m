function [Wx, Wy, varargout] = sinusoidal_wind(x, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		sinusoidal_wind
%
% PURPOSE:		Calculate a wind field at a set of locations in space
%            
% SYNTAX:		[Wx, Wy] = sinusoidal_wind(x)
%				[Wx, Wy] = sinusoidal_wind(x, V, lambda)
%				[Wx, Wy, dWx_dx, dWx_dy] = sinusoidal_wind(x, V, lambda)
%
% INPUTS:		x		- x coordinate (along wind direction)
%				V		- Wind speed
%				lambda	- Wind field wavelength
%
% OUTPUTS:		Wx		- x wind velocity
%				Wy		- y wind velocity
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		December 2008
%
% MODIFIED:		December 2008
%
% See also:		sinusoidal_wind3, linear_wind
%--------------------------------------------------------------------------
if nargin == 1
	V = 5;
	lambda = 10;
else
	V = varargin{1};
	lambda = varargin{2};
end


theta = pi/4*cos(2*pi*x/lambda);

Wx = -V*cos(theta);
Wy = -V*sin(theta);

if nargout >2
	dtheta_dx = -pi*pi/(2*lambda)*sin(2*pi*x/lambda);
	varargout{1} =-dtheta_dx.*Wy;
	varargout{2} = dtheta_dx.*Wx;
end
