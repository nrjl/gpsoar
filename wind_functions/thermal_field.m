function [W, Jw] = thermal_field(X, thermal_params)
%--------------------------------------------------------------------------
%
% FUNCTION:		thermal_field
%
% PURPOSE:		Calculate 3D wind field and local gradients for a field of 
%				toroidal thermals
%
% SYNTAX:		[W, Jw] = torus_thermal(pos, thermal_params)
%
% INPUTS:		pos		- [3×n] matrix of coordinates
%				params  - [k×6] each row containing thermal parameters:
%							[w_core, R, k, centre(x,y,z)]
%
% OUTPUTS:		W		- [3×n] matrix of x,y,z components of wind
%							Note that wind is zero outside thermal region
%				Jw		- [3×3×n] matrix of local gradients matrices
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		January 2009
%
% MODIFIED:     January 2009
%
% See also:		torus_thermal3
%--------------------------------------------------------------------------
n_thermals = size(thermal_params,1);
W = zeros(size(X));
Jw = zeros(3,3,size(X,2));

for i = 1:n_thermals
	R = thermal_params(i,2); k = thermal_params(i,3);
	inrange = (X > (thermal_params(4:6)'- [2*R;2*R;k*R])*ones(1,size(X,2)));
	inrange = inrange & (X < (thermal_params(4:6)'+ [2*R;2*R;k*R])*ones(1,size(X,2)));
	
	if any(all(inrange))
		[W1, Jw1] = torus_thermal3(X, thermal_params(i,:));
		W = W + W1;
		Jw = Jw + Jw1;
	end
end