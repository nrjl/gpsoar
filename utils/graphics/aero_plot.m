function h_aero = aero_plot(plane_aero)
%--------------------------------------------------------------------------
%
% FUNCTION:		aero_plot
%
% PURPOSE:		plot the aerodynamic surfaces of an aircraft
%
% SYNTAX:		aero_plot(plane_aero, X)
%
% INPUTS:		plane_aero	- structure containing definition of
%					aerodynamic surfaces of an aircraft (See SBXC_def)
%				X			- state vector
%
% OUTPUTS:		(plots to current figure)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2007
%
% MODIFIED:     October 2007
%
% See also:		strip_plot, plane_movie, SBXC_def
%--------------------------------------------------------------------------

n_surf = length(plane_aero);

h_aero = hggroup;

for i = 1:n_surf
	qc = plane_aero(i).qc_surf;
	ai = plane_aero(i).base_alfa;
	mirror = plane_aero(i).mirror;
	
	n_poly = size(qc,2)-1;
	
	for j = 1:n_poly
		xp = [qc(1,j)  + 0.25*qc(4,j), qc(1,j+1)+ 0.25*qc(4,j+1), ...
			  qc(1,j+1)- 0.75*qc(4,j+1), qc(1,j)  - 0.75*qc(4,j)];
			   
		yp = [qc(2,j), qc(2,j+1), qc(2,j+1), qc(2,j)];
		zp = [qc(3,j), qc(3,j+1), qc(3,j+1), qc(3,j)];
		
		poly = [xp; yp; zp];
		h_poly = fill3(poly(1,:), poly(2,:), poly(3,:), 'g');
		set(h_poly, 'Parent', h_aero);
		
		if mirror
			poly = [xp; -yp; zp];
			h_poly = fill3(poly(1,:), poly(2,:), poly(3,:), 'r');
			set(h_poly, 'Parent', h_aero);
		end
	end
	
end