function u_log = log_profile(z, u_inf, delta, z0)
%--------------------------------------------------------------------------
%
% FUNCTION:		log_profile
%
% PURPOSE:		Calculate horizontal wind value in a logarithmic profile
%
% SYNTAX:		u_log = log_profile(z, u_inf, delta, z0)
%
% INPUTS:		z		- altitude vector above layer start
%				u_inf	- free stream velocity (above layer)
%				delta	- layer thickness
%				z0		- Roughness length
%
% OUTPUTS:		u		- horizontal wind velocities
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		June 2007
%
% MODIFIED:     June 2007
%
% See also:		profile_flythrough, log_profile, wind_field
%--------------------------------------------------------------------------
% Logarithmic

zd = (z<delta).*z.*(z>=z0)+z0*(z<z0 | z>=delta);

u_log = u_inf*(log(zd/z0)/log(delta/z0)) + (z>=delta)*u_inf;