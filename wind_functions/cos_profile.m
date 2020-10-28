function u = cos_profile(z, zl, zu, vl, vu)
%--------------------------------------------------------------------------
%
% FUNCTION:		cos_profile
%
% PURPOSE:		Calculate horizontal wind value in a cosine profile
%
% SYNTAX:		u = cos_profile(z, zl, zu, vl, vu)
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
% See also:		log_profile, wind_field
%--------------------------------------------------------------------------

z_cos = (z > zl) & (z < zu);
u = (z <= zl)*vl + z_cos.*((vu-vl)/2*(1 - cos((z-zl)/(zu-zl)*pi)) + vl) + (z >= zu)*vu;