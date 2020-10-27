function [SBXC_aero, SBXC_param] = SBXC_def(varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		SBXC_def
%
% PURPOSE:		define aerodynamic surfaces for the SBXC glider
%               
% SYNTAX:		[SBXC_aero, m, Ixx, Iyy, Izz, Ixz] = SBXC_def()
%				... = SBXC_def(x_cg)
%				... = SBXC_def(x_cg, y_cg)
%				... = SBXC_def(x_cg, y_cg, z_cg)
%
% INPUTS:		x_cg, y_cg, z_cg - CG position relative to wing leading
%					edge, in standard aircraft axes (forward, right wing,
%					down)
%
% OUTPUTS:		SBXC_aero	- Structure contiaing following fields for each
%						surface:
%					name		- string containing surface name
%					qc_surf		- [4×n] array containg quarter chord 
%								positions and chord length in each column
%					base_alfa	- base angle of attack (rad)
%					mirror		- logical for mirrored surfaces (through xz
%									plane)
%					profile		- {1×2} cell array containing coefficient
%								table and Reynolds number (see plot_polars)
%
%				m			- Total aircraft mass (kg)
%				Ixx, Iyy, Izz, Ixz - Moments of inertia
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		September 2007
%
% MODIFIED:     October 2007
%
% See also:		strip_method, strip_forces, aero_plot, plot_polars
%--------------------------------------------------------------------------

%% Input check
root_chord = 290;

switch nargin
	case 0
		x_cg = -0.4*root_chord;
		y_cg = 0;
		z_cg = 0;
	case 1
		x_cg = varargin{1};
		y_cg = 0;
		z_cg = 0;
	case  2
		x_cg = varargin{1};
		y_cg = varargin{2};
		z_cg = 0;
	case  3
		x_cg = varargin{1};
		y_cg = varargin{2};
		z_cg = varargin{3};
	otherwise
		disp('Too many input arguments, reverting to default')
		x_cg = -0.4*root_chord;
		y_cg = 0;
		z_cg = 0;
end



%% AIRCRAFT MASS PROPERTIES AND FUSELAGE FINENESS RATIO
SBXC_param.m = 5.443;
SBXC_param.Ixx = 1374*(4.3/15)^2*5.443/363;
SBXC_param.Iyy = 869*(4.3/15)^2*5.443/363;
SBXC_param.Izz = 2214*(4.3/15)^2*5.443/363;
SBXC_param.Ixz = 67*(4.3/15)^2*5.443/363;
SBXC_param.fuse_l = 25*3*25.4/1000;			% Fuselage length
SBXC_param.fuse_d = 24/3*25.4/1000;			% Fuselage diameter
SBXC_param.fuse_cg = 0.25;					% Fuselage cg location/fuse_l
SBXC_param.Cd0 = 0.05;						% Estimated profile drag

%% AIRCRAFT GEOMETRY
% Define quarter chords
load s2048_control.mat
s2048_control = {polars, Re, def};

load NACA0010-54_control.mat
NACA0010_control = {polars, Re, def};

% Main Wing (right)
wing_incidence = 3.58;

% Inboard section (with flaps)
wing_flaps = [-290/4, -65-225/4;...		% x coord of quarter chord
				100,     1356.8;...		% y coord of quarter chord
				0,     -90.1;...		% z coord of quarter chord
			  290,      225];			% chord length

wing_flaps = wing_flaps - [x_cg; y_cg; z_cg; 0]*ones(1, size(wing_flaps, 2));
qc_wing_flaps = segmenter(wing_flaps, 6)/1000;	% In body coordinates
control_flaps = 4*ones(1, size(qc_wing_flaps, 2)-1);
reverse_flaps = 0*ones(1, size(qc_wing_flaps, 2)-1);

% Mid section (with ailerons)
wing_ailrn = [-65-225/4, -144-143.5/4;...
			1356.8,    2100.3;...
			-90.1,    -227.1;...
			225,      143.5];

wing_ailrn = wing_ailrn - [x_cg; y_cg; z_cg; 0]*ones(1, size(wing_ailrn, 2));
qc_wing_ailrn = segmenter(wing_ailrn, 5)/1000;	% In body coordinates
control_ailrn = 2*ones(1, size(qc_wing_ailrn, 2)-1);
reverse_ailrn = -1*ones(1, size(qc_wing_ailrn, 2)-1);

% Tip section (no control)
wing_tip = [-144-143.5/4, -278.5-11.5/4;...
           2100.3,       2246.9;...
           -227.1,       -253.5;...
           143.5,         11.5];

wing_tip = wing_tip - [x_cg; y_cg; z_cg; 0]*ones(1, size(wing_tip, 2));
qc_wing_tip = segmenter(wing_tip, 1)/1000;		% In body coordinates
control_tip = 0*ones(1, size(qc_wing_tip, 2)-1);
reverse_tip = 0*ones(1, size(qc_wing_tip, 2)-1);

SBXC_aero(1).name		= 'Main wing';
SBXC_aero(1).qc_surf	= [qc_wing_flaps, qc_wing_ailrn(:,2:end), qc_wing_tip(:,2:end)];
SBXC_aero(1).base_alfa	= wing_incidence*ones(1, size(SBXC_aero(1).qc_surf, 2)-1);
SBXC_aero(1).mirror		= 1;
SBXC_aero(1).profile	= s2048_control;
SBXC_aero(1).control	= [control_flaps, control_ailrn, control_tip];
SBXC_aero(1).reverse	= [reverse_flaps, reverse_ailrn, reverse_tip];
SBXC_aero(1).c_lim		= [[7.79; -7.79], [4.55; -7.6], [21.28; -21.28], [90; -5]];

% Stabilator
r_stab = [-(1170+125/4), -(1194.5+100/4);...
          0, 400;...
          -70, -70;...
          125, 100];         
 
r_stab = r_stab - [x_cg; y_cg; z_cg; 0]*ones(1, size(r_stab, 2));
qc_stab = segmenter(r_stab, 3)/1000;

SBXC_aero(2).name		= 'Horizontal stabilator (elevators)';
SBXC_aero(2).qc_surf	= qc_stab;
SBXC_aero(2).base_alfa	= zeros(1, size(qc_stab, 2)-1);
SBXC_aero(2).mirror		= 1;
SBXC_aero(2).profile	= NACA0010_control;
SBXC_aero(2).control	= -1*ones(1, size(qc_stab,2)-1);	% This denotes all-moving surface
SBXC_aero(2).reverse	= zeros(1, size(qc_stab,2)-1);
SBXC_aero(2).c_lim		= SBXC_aero(1).c_lim;

% Fin
fin = [-(1104.5+255/4), -(1318.1 + 165/4);...
        0, 0;...
        10, -360;...
        255, 165];

fin = fin - [x_cg; y_cg; z_cg; 0]*ones(1, size(fin, 2));
fin_sections = 4;
qc_fin = segmenter(fin, fin_sections)/1000;

SBXC_aero(3).name		= 'Fin (rudder)';
SBXC_aero(3).qc_surf	= qc_fin;
SBXC_aero(3).base_alfa	= zeros(1, size(qc_fin, 2)-1);
SBXC_aero(3).mirror		= 0;
SBXC_aero(3).profile	= NACA0010_control;
SBXC_aero(3).control	= 3*ones(1, size(qc_fin,2)-1);
SBXC_aero(3).reverse	= zeros(1, size(qc_fin,2)-1);
SBXC_aero(3).c_lim		= SBXC_aero(1).c_lim;