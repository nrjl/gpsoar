function strip_plot(X, U)
%--------------------------------------------------------------------------
%
% FUNCTION:		strip_plot ****** DEAD CODE **********
%				Replaced by strip_method with plot flag
%
% PURPOSE:		calculate and plot aircraft body force and moments using a 
%				simple strip method
%               
% SYNTAX:		strip_plot(X, U)
%
% INPUTS:		X		- state vector
%				U		- control vector
%
% OUTPUTS:		poa_all	- Points of application in body axes
%				F_all	- All forces in body axes
%				M_all	- All moments in body axes
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2007
%
% MODIFIED:     October 2007
%
% See also:		strip_forces, plane_movie
%--------------------------------------------------------------------------
global PLANE_AERO

% stime = cputime;
l_conv = eye(4);
l_conv(2,2) = -1;

n_surf = length(PLANE_AERO);

k = 1;
n_poa = 0;
cn = 1;

for i = 1:n_surf
	n_poa = n_poa + size(PLANE_AERO(i).qc_surf, 2)*(1 + PLANE_AERO(i).mirror);
end

poa_all = zeros(3, n_poa);
F_all = poa_all;
M_all = poa_all;

for i = 1:n_surf
	qc_surf = PLANE_AERO(i).qc_surf;
	ai = PLANE_AERO(i).base_alfa;
	mirror = PLANE_AERO(i).mirror;
	profile = i;
	control = PLANE_AERO(i).control;
	
	if control; def = U(control); else def = 0; end
	
	fn = cn + size(PLANE_AERO(i).qc_surf, 2)-2;

	[poa_all(:, cn:fn), F_all(:, cn:fn), M_all(:, cn:fn)] = ...
		strip_forces(qc_surf, X, ai, @wing_lookup, profile{1}, profile{2}, def);
	
	k = k+1;
	cn = fn+1;
	
	if mirror
		qc_surf = l_conv*qc_surf;
		
		fn = cn + size(PLANE_AERO(i).qc_surf, 2)-2;
		[poa_all(:, cn:fn), F_all(:, cn:fn), M_all(:, cn:fn)] = ...
			strip_forces(qc_surf, X, ai, @wing_lookup, profile{1}, ...
			profile{2}, def);
		
		k = k+1;
		cn = fn+1;
	end	
end

% max_F = max(sqrt(F_all(1,:).^2 + F_all(2,:).^2 + F_all(3,:).^2));
% F_scale = .290/max_F;
F_scale = 0.1;

Ceb = calc_Ceb(X(7:9));

F_all = Ceb*F_all*F_scale;
F_end = poa_all + F_all;

aw = sqrt(sum(F_all.^2, 1))*1.5;	% Arrow width

% for i = 1:size(F_end, 2)
% 	plot3([poa_all(1,i), F_end(1,i)], [poa_all(2,i), F_end(2,i)], ...
% 		[poa_all(3,i), F_end(3,i)], 'b--');
% end

arrow3(poa_all', F_end', 'b-', aw, 3*aw, 0.2);