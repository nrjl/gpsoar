function axis_lim = plane_plot(pts, x, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		plane_plot
%
% PURPOSE:		Plot the geometry for an aeroplane model
%               
% SYNTAX:		axis_lim = plane_plot(pts, x, dot_flag)
%
% INPUTS:		pts     - aeroplane points (structure with a 3xn array for
%                         each component (normally wing, tail plane, fin)
%				x		- state vector
%						  [u, v, w, p, q, r, psi, theta, phi, x, y, z]
%               dot_flag- flag to determine whether dots should be drawn at
%                         each point
%
% OUTPUTS:		axis_lim- axis limits
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		March 2007
%
% MODIFIED:     June 2007
%
% See also:		sbxc_full, plot_polars
%--------------------------------------------------------------------------

hold on;
n_components = size(pts, 2);

Ceb = calc_Ceb(x(7:9));
	
axis_lim = [0; 0; 0];
% x(12) = -x(12);

for i = 1:n_components
	a = pts{i};
	a1 = Ceb*a/1000 + x(10:12)*ones(1, size(a, 2));
	a2 = Ceb*([1,0,0;0,-1,0;0,0,1]*a)/1000 + x(10:12)*ones(1, size(a, 2));

    fill3(a1(1,:), a1(2,:), a1(3,:), 'r')    
    fill3(a2(1,:), a2(2,:), a2(3,:), 'r')
    
    max_a = max(abs(a), [], 2);
    axis_lim = max(max_a, axis_lim);

    if nargin > 2
		if varargin{1}
			plot3(a(1,:), a(2,:), a(3,:), 'k.', 'MarkerSize', 5)
			plot3(a(1,:), -a(2,:), a(3,:), 'k.', 'MarkerSize', 5)
		end
    end
end

% [xs, ys, zs] = sphere(20);
% xs = (xs*0.02 + x(10)*ones(size(xs)));
% ys = (ys*0.02 + x(10)*ones(size(ys)));
% zs = (zs*0.02 + x(10)*ones(size(zs)));
% surf(xs, ys, zs);

% plot3([0 200], [0 0], [0 0], 'b-'); 
% plot3([0 0], [0 200], [0 0], 'b-')
% plot3([0 0], [0 0], [0 200], 'b-')
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
axis equal;
axis tight;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(50, 18);
% xy_lim = max(axis_lim(1:2));
% axis([-xy_lim, xy_lim, -xy_lim, xy_lim, -axis_lim(3), axis_lim(3)])