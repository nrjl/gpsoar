function [W, varargout] = torus_thermal3(pos, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		torus_thermal3
%
% PURPOSE:		Calculate 3D wind field and local gradients for toroidal 
%				thermal model
%
% SYNTAX:		[W] = torus_thermal(pos, w_core, R, k, centre)
%				[W] = torus_thermal(pos, params)
%				[W, Jw] = torus_thermal(pos, w_core, R, k, centre)
%
% INPUTS:		pos		- [3×n] matrix of coordinates
%				params  - [1×6] vector containing [w_core, R, k, centre]
%				w_core	- core vertical velocity
%				R		- major radius
%				k		- elliptical factor k = (major axis)/(minor axis)
%				centre	- [3×1] vector of thermal centre
%
% OUTPUTS:		W		- [3×n] matrix of x,y,z components of wind
%							Note that wind is zero outside thermal region
%				Jw		- [3×3×n] matrix of local gradients matrices
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2007
%
% MODIFIED:     November 2007
%
% See also:		torus_thermal
%--------------------------------------------------------------------------
if nargin == 2
	w_core	= varargin{1}(1);
	R		= varargin{1}(2);
	k		= varargin{1}(3);
	centre	= varargin{1}(4:6);
elseif nargin == 5
	w_core	= varargin{1};
	R		= varargin{2};
	k		= varargin{3};
	centre	= varargin{4};
end
	

W = torus_thermal(pos, w_core, R, k, centre(:));
dx = 0.1;

if nargout == 2	
	Jw = zeros(3,3,size(pos,2));
	for i = 1:3
		pos2 = pos;
		pos2(i,:) = pos(i,:)+dx;
		Jw(i,:,:) = ((torus_thermal(pos2, w_core, R, k, centre(:)) - W)/dx);
	end
	varargout{1} = Jw;
end