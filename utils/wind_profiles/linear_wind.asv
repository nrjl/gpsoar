function [Wx, Wy] = linear_wind(x, y, Wx0, Wy0, varargin)
%--------------------------------------------------------------------------
% FUNCTION:		linear_wind
%
% PURPOSE:		Calculate wind in 2D form from linear derivatives
%
% SYNTAX:		[Wx, Wy] = linear_wind(x, y, Wx0, Wy0, Jw)
%				[Wx, Wy] = linear_wind(x, y, Wx0, Wy0, dWx_dx, dWx_dy,
%								dWy_dx, dWy_dy)
%
% INPUTS:		x	- [1×n] vector of x positions
%				y	- [1×n] vector of y positions
%				Wx0	- [1×1] Wind in x direction at origin (0,0)
%				Wy0	- [1×1] Wind in y direction at origin (0,0)
%				Jw	- [2×2] Matrix of wind derivatives: [dWx_dx, dWx_dy
%														 dWy_dx, dWy_dy]
%
% OUTPUTS:		Wx	- x wind at specified locations
%				Wy	- y wind at specified locations
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2008
%
% MODIFIED:     January 2009
%
% See also:		wind_field, sinusoidal_wind
%--------------------------------------------------------------------------
if nargin > 5
	dWx_dx = varargin{1};
	dWx_dy = varargin{2};
	dWy_dx = varargin{3};
	dWy_dy = varargin{4};
else
	Jw = varargin{1};
	dWx_dx = Jw(1,1);
	dWx_dy = Jw(1,2);
	dWy_dx = Jw(2,1);
	dWy_dy = Jw(2,2);
end	

Wx = Wx0 + dWx_dx*x + dWx_dy*y;
Wy = Wy0 + dWy_dx*x + dWy_dy*y;