function varargout = torus_thermal(varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		torus_thermal
%
% PURPOSE:		calculate wind field for toroidal thermal model
%
% SYNTAX:		V_wg_e = torus_thermal(loc, w_core, R, k, centre)
%				[U,V,W] = torus_thermal(X, Y, Z, w_core, R, k, centre)
%
% INPUTS:		loc		- [3×n] matrix of coordinates
%				w_core	- core vertical velocity
%				R		- major radius
%				k		- elliptical factor k = (major axis)/(minor axis)
%				centre	- [3×1] vector of thermal centre
%
% OUTPUTS:		V_wg_e	- [3×n] matrix of x,y,z components of wind
%							Note that wind is zero outside thermal region
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2007
%
% MODIFIED:     November 2007
%
% See also:		wind_field, ellipse_field
%--------------------------------------------------------------------------

if nargin == 5
	loc = varargin{1};
	X = loc(1,:);
	Y = loc(2,:);
	Z = loc(3,:);
	w_core = varargin{2};
	R = varargin{3};
	k = varargin{4};
	centre = varargin{5};
elseif nargin == 7
	X = varargin{1};
	Y = varargin{2};
	Z = varargin{3};
	w_core = varargin{4};
	R = varargin{5};
	k = varargin{6};
	centre = varargin{7};
else
	error('thermal_model:torus_thermal:nargin',...
		'Incorrect number of input arguments')
end

X = X-centre(1);
Y = Y-centre(2);
Z = Z-centre(3);

d = sqrt(X.^2 + Y.^2);

u2 = atan2(Y,X);
% t2 = atan2(Z,(d-R)*k^2).*(X|Y) + pi*~(X|Y);

zerod = (d == 0);
% bigd = (d < 2*R) & (abs(Z) < (k*R*2));

% rall = 1/k*sqrt((d - R).^2*k^2  + Z.^2);
rall = sqrt(X.^2 + Y.^2 + (2*Z/k).^2);
rcut = rall <= 2*R;

d = d+zerod;
dR0 = (d-R == 0);
dR = d - R + dR0;

height_adjust = 0.5.*(cos(pi*Z/(k*R)) + 1);

W = -w_core*R/pi*sin(pi*d/R)./d;
W = W.*(~zerod & rcut) - (zerod & rcut)*w_core;
W = W.*height_adjust;
Wd = (-W.*Z./(dR*k^2)).*~dR0 - w_core*Z./(R*k^2).*height_adjust.*dR0.*rcut;
U = Wd.*cos(u2);
V = Wd.*sin(u2);

if nargout == 1
	varargout{1} = [U; V; W];
elseif nargout == 3
	varargout{1} = U;
	varargout{2} = V;
	varargout{3} = W;
elseif nargout ~= 0
	error('thermal_model:torus_thermal:nargout',...
		'Incorrect number of output arguments')
end



% a = -w_core*R/pi*sin(pi*d/R)./(d.*cos(t2));
% a = a.*(~zerod & rcut) + (zerod & rcut)*w_core;
% a = .5*a.*(cos(pi*Z/(k*R)) + 1);
% 
% U = -a.*sin(t2).*cos(u2);
% V = -a.*sin(t2).*sin(u2);
% W = a.*cos(t2);