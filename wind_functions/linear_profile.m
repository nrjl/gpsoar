function u = linear_profile(z, zl, zu, vl, vu)
%--------------------------------------------------------------------------
%
% FUNCTION:		linear_profile
%
% PURPOSE:		Calculate horizontal wind value in a linear profile
%
% SYNTAX:		u = linear_profile(z, zl, zu, vl, vu)
%
% INPUTS:		z	- altitude vector
%				zl	- lower altitude limit of profile
%				zu	- upper altitude limit of profile
%				vl	- lower velocity limit of profile
%				vu	- upper velocity limit of profile
%
% OUTPUTS:		u	- horizontal wind velocities
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2008
%
% MODIFIED:
%
% See also:		log_profile, cos_profile, wind_field
%--------------------------------------------------------------------------

z_lin = (z > zl) & (z < zu);
u = (z <= zl)*vl + z_lin.*((z-zl)*(vu-vl)/(zu-zl) + vl) + (z >= zu)*vu;