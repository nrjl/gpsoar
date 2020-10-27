function varargout = plane_properties(aero, param)
%--------------------------------------------------------------------------
%
% FUNCTION:		plane_properties
%
% PURPOSE:		Display the basic aerodynamic parameters of an aeroplane 
%               model
%
% SYNTAX:		plane_properties(aero, param)		Display properties
%				S = plane_properties(aero, param)	Return reference area
%				[S, AR] = plane_properties(aero, param)
%
% INPUTS:		aero	- aerodynamic definition structure (from SBXC_def)
%				param	- aircraft parameters
%
% OUTPUTS:		(display to command window)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2007
%
% MODIFIED:     November 2007
%
% See also:		SBXC_def
%--------------------------------------------------------------------------
fprintf('\n----- Aerodynamic Components -----\n')

n_components = length(aero);
c_area = zeros(1, n_components);
b = c_area;

for i = 1:n_components
	name = aero(i).name;
	qc = aero(i).qc_surf;
	c_area(i) = sqrt(sum(diff(qc(1:3,:),1,2).^2))*...
		(qc(4,1:end-1) + diff(qc(4,:))*0.5)';
	c_area(i) = c_area(i)*(1 + aero(i).mirror);
	b(i) = sqrt(sum((qc(2:3,end) - qc(2:3,1)).^2))*(1 + aero(i).mirror);
	fprintf(1, 'Component %i: %s\n\tSpan = %0.5g m\n', i, name, b(i));
	fprintf(1, '\tArea = %0.5g m^2\n\tAspect ratio = %0.5g\n', c_area(i), b(i)*b(i)/c_area(i));
end
fprintf(1, 'Total aerodynamic area = %0.4g m^2\n', sum(c_area));

fprintf(1, '\n----- Fuselage -----\n')
fprintf(1, 'Length = %0.5g m\nDiameter = %0.5g m\nCd0 = %0.5g\n', ...
	param.fuse_l, param.fuse_d, param.Cd0)


fprintf(1, '\n----- Mass Components -----\n')
fprintf(1, 'Total vehicle mass = %0.5g\n', param.m);
fprintf(1, ['Moments of Inertia (kg.m^2)\n\tIxx = %0.5g\tIyy = %0.5g\n\tIzz = ', ...
	'%0.5g\tIxz = %0.5g\n'], param.Ixx, param.Iyy, param.Izz, param.Ixz)

if nargout == 1
	varargout{1} = c_area(1);
elseif nargout == 2
	varargout{1} = c_area(1);
	varargout{2} =  b(1)*b(1)/c_area(1);
end