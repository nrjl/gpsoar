function K = projected_square_expt(x1, t1, y1, x2, t2, loghyper, xdist, tdist, pdist)
% loghyper = log([lx, lt, lp, sigma_f (, sigma_n)]
if nargin < 9		%nargin = 3 or 5 -> calculate sqdist
	xdist = square_dist(x1, x2);
	tdist = square_dist(t1, t2);
	pdist = projected_point_distance(x1, t1, y1, x2, t2);
end

K = exp(2*loghyper(4)-0.5/exp(2*loghyper(1))*xdist ...
					- 0.5/exp(2*loghyper(2))*tdist ...
					- 0.5/exp(2*loghyper(3))*pdist);