function h_plane = SBXC_handle(X, plane_aero, varargin)


hold on;

h_fuse = SBXC_fuse( );
h_aero = aero_plot(plane_aero);

h_plane = hgtransform;
set(h_fuse, 'Parent', h_plane);
set(h_aero, 'Parent', h_plane);

Ceb = calc_Ceb(X(7:9));
Cr = eye(4);
Cr(1:3, 1:3) = Ceb;			% Rotation transformation

Ct = eye(4);
Ct(1:3,4) = X(10:12);		% Translate to actual position

set(h_plane, 'Matrix', Ct*Cr)

if nargin == 3
	set(gca,'YDir','reverse');
	set(gca,'ZDir','reverse');
	axis equal;
	grid on;
	xlabel('X'); ylabel('Y'); zlabel('Z');
	view(50, 18);
end