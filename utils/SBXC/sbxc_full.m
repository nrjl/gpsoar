function full_plane = sbxc_full(varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		sbxc_full
%
% PURPOSE:		Calculate the geometry for the SBXC Remote Control Glider
%               
% SYNTAX:		full_plane = sbxc_full()		standard cg position
%				full_plane = sbxc_full(x_cg)
%				full_plane = sbxc_full(x_cg, y_cg)
%				full_plane = sbxc_full(x_cg, y_cg, z_cg)
%
% INPUTS:		x_cg, y_cg, z_cg - CG position relative to wing leading
%					edge, in standard aircraft axes (forward, right wing,
%					down)
%
% OUTPUTS:		full_plane	- {1 × 3} cell array for wing, tail, and fin
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		March 2007
%
% MODIFIED:     June 2007
%
% See also:		plane_plot, plot_polars
%--------------------------------------------------------------------------

%% Introduction
% Standard aircraft axes definition; x forward, y starboard, z down
aircraft_name = 'SB/XC Glider';

%% Primary (inner) wing section parameters
root_chord = 290;
inner_span = 1260;

switch nargin
	case 0
		x_cg = -0.5*root_chord;
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

section_chord = 225;
inner_n = 20;

root_dihedral = 4.1; % Degrees

y_inner = linspace(0, inner_span, inner_n);
c_inner = linspace(root_chord, section_chord, inner_n);

inner_area = trapz(y_inner, c_inner);

inner_pts = [c_inner; y_inner; zeros(1, inner_n)];
inner_pts = rotate_x(-root_dihedral, 'deg')*inner_pts;

%% Outer span parameters
outer_span = 905;
tip_dihedral = 6.3; % degrees

y_tip = [0, 2, 4, 6, 8, 10:8:266, 282:16:(282+12*16)]; % From tip
y_tip = reverse(outer_span - y_tip); % From root

c_tip = [0 13.5 22 28 31.5 37 53 67 77 86.5 94.5 102 108 ...
114 119 123 127.5 131 134 137.5 140.5 144 146 149 151 154 ...
155.5 158 160 161.5 163 164.5 166 167.5 169 170 171 172 174.5 ...
176.5 178 180 181.5 183 184 185.5 186.5 188 189.5 191 192.5];
c_tip = reverse(c_tip);

outer_n_pts = 10;
y_mid = linspace(0, min(y_tip), outer_n_pts);
c_mid = linspace(section_chord, max(c_tip), outer_n_pts);

y_outer = [y_mid, y_tip];
c_outer = [c_mid, c_tip];
outer_area = trapz(y_outer, c_outer);

outer_pts = [c_outer; y_outer; zeros(1,length(c_outer))];
outer_pts = rotate_x(-(tip_dihedral+root_dihedral), 'deg')*outer_pts;

offset = inner_pts(:,inner_n);
offset(1,1) = 0; % Do not want to add chord offset

outer_pts = outer_pts + offset*ones(1,length(c_outer));

%% Full wing points

wing_ac = 0.25;     % Assume quarter-chord A/C
wing_incidence = 3.48; % Angle of incidence (deg)
x_cg = x_cg + wing_ac*root_chord;


% ----- Flattened wing points ----- %
% y_full = [y_inner, y_outer];
% c_full = [c_inner, c_outer];
wing_area = (inner_area + outer_area)/1e6;

% ----- Trailing edge definition ----- %
mid_trail = [0; inner_pts(2,inner_n); inner_pts(3,inner_n)];
root_trail = [0; 0; 0;];

wing_pts = [inner_pts, outer_pts, mid_trail, root_trail];

wing_pts(1,:) = -(1 - wing_ac)*root_chord + wing_pts(1,:);

wing_pts = rotate_y(wing_incidence, 'deg')*wing_pts;

wing_pts(1,:) = -x_cg + wing_pts(1,:);
wing_pts(2,:) = -y_cg + wing_pts(2,:);
wing_pts(3,:) = -z_cg + wing_pts(3,:);


%% Tailplane (elevator) Definition
% All moving tailplane/elevator

tail_xle = -x_cg - 1170 - wing_ac*cos(wing_incidence)*root_chord;
tail_yle = -y_cg + 0;
tail_zle = -z_cg - 70;

tail_root_chord = 124.5;
tail_tip_chord = 100;
tail_halfspan = 400;

tail_n_pts = 5;
tail_y = [linspace(0, tail_halfspan, tail_n_pts), tail_halfspan, 0];
tail_c = [linspace(tail_root_chord, tail_tip_chord, tail_n_pts), 0, 0];

tail_pts = [(tail_xle + tail_c); (tail_yle + tail_y); ...
    tail_zle*ones(1,length(tail_y))];

tail_area = trapz(tail_y, tail_c)/1e6;

%% Fin Definition

fin_sweep = 60; % deg
fin_root_chord = 255;
fin_tip_chord = 165;
fin_height = 370;

fin_xle = -x_cg - 980 - wing_ac*cos(wing_incidence)*root_chord;
fin_yle = -y_cg + 0;
fin_zle = -z_cg + 10;

fin_a = fin_height/tand(fin_sweep);
fin_pts = [fin_xle-[0, fin_a, fin_a + fin_tip_chord, fin_root_chord]; ...
           fin_yle+[0, 0, 0, 0]; ...
           fin_zle-[0, fin_height, fin_height, 0] ];
       
fin_area = (fin_root_chord+fin_tip_chord)/2*fin_height/1e6;
       

%% Fuselage
fuse_length = 25*3*25.4;
r_fuse = [0, 4, 7, 9, 11, 11.5, 12, 12, 11.5, 11, 10, 9, 7, 6, 5.5, 5, ...
	4.5, 4, 3.5, 3, 3, 2.5, 2, 1.5, 1, 0.5]/3*25.4;


%% Output
% ----- Command window output ----- %
fprintf(1, '\n ----- %s ----- \n', aircraft_name)

fprintf(1, 'Main wing area\t= %0.4g m^2\n', wing_area*2)
fprintf(1, 'Tailplane area\t= %0.5g m^2\n', tail_area*2)
fprintf(1, 'Tail fin area\t= %0.5g m^2\n', fin_area)
fprintf(1, 'Reference area\t= %0.5g m^2\n', 2*(wing_area+tail_area))
fprintf(1, 'Reference chord\t= %0.5g m\n', section_chord/1e3)
fprintf(1, 'Reference span\t= %0.5g m\n', 2*max(wing_pts(2,:))/1e3)

% ----- Full Aircraft Structure ----- %
% aero = struct('surface', {wing_pts, tail_pts, fin_pts}, 'mirror', {1, 1, 0});
% fuse = struct('radii', {r_fuse}, 'length', {fuse_length}, 'tip_loc', {tip_loc});
% nplane = struct('aero', {aero}, 'fuse', {fuse});

% ----- Figures ----- %
full_plane = {wing_pts, tail_pts, fin_pts};
% x = zeros(1, 12);
% 
% figure(4); clf;
% hold on;
% 
% subplot(2,2,1);
% axis_lim = plane_plot(full_plane, x);
% view(33,22);
% axis off;
% 
% subplot(2,2,2);
% plane_plot(full_plane, x);
% view(90,90);
% grid on;
% 
% subplot(2,2,3);
% plane_plot(full_plane, x);
% view(0,0);
% grid on;
% 
% subplot(2,2,4);
% plane_plot(full_plane, x);
% view(90,0);
% grid on;
