function [u, varargout] = pohlhausen(z, z0, u_inf, delta, lambda)
%--------------------------------------------------------------------------
%
% FUNCTION:		pohlhausen
%
% PURPOSE:		Calculate horizontal wind value in a pohlhausen's quartic
%				profile
%
% SYNTAX:		u = pohlhausen(z, z0, u_inf, delta, lambda)
%
% INPUTS:		z		- target altitude points
%				z0		- layer start altitude (lower limit)
%				u_inf	- free stream velocity (above layer)
%				delta	- layer thickness
%				lambda	- Pohlhausen's quartic parameter (normally 1)
%
% OUTPUTS:		u		- horizontal wind velocities
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		June 2007
%
% MODIFIED:     March 2008
%
% See also:		profile_flythrough, log_profile, wind_field
%--------------------------------------------------------------------------
z = z-z0;
zd = (z<=delta).*z.*(z>=0);

% Pohlhausen's quartic
eta = zd/delta;
u = u_inf*(2*eta - 2*eta.^3 + eta.^4 + lambda/6*eta.*(1-eta).^3 + (z>delta));

if nargout > 1
	varargout{1} = u_inf*(2/delta - 6*zd.^2/delta^3 + 4*zd.^3/delta^4 + ...
		lambda/(6*delta)*(1 - eta).^2.*(1 + 2*eta)).*(z>0 & z<delta);
end