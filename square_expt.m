function Kxt = square_expt(x1, t1, x2, t2, loghyper, d_xx, d_tt)
% loghyper = log([l, lt, sigma_f, sigma_n]);

n_dim = size(x1, 2)*(numel(loghyper)>4) + (numel(loghyper) <= 4);

if nargin < 7
	d_tt = square_dist(t1, t2);
	if n_dim > 1
		d_xx = square_dist(x1, x2, exp(loghyper(1:n_dim)));
	else
		d_xx = square_dist(x1, x2);
	end
	
end

Kxx = square_exp(x1, x2, loghyper([1:n_dim,n_dim+2]), d_xx);		% Training spatial covariances
Ktt = square_exp(t1, t2, [loghyper(n_dim+1), 0], d_tt);		% Training temporal covariances

Kxt = Kxx.*Ktt;