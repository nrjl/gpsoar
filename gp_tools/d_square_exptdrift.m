function C = d_square_exptdrift(x, t, index, loghyper)
% Different from Rasmussen!!
% loghyper: log([length_scale, lag_scale, sigma_f, sigma_n])
ll = loghyper(1);
llt = loghyper(2);
lsf = loghyper(3);

if index == 1		% dK/dl
	sq_dist_xx = square_dist(x, x);
	sq_dist_tt = square_dist(t, t);
	C = exp(2*lsf - 2*ll - sq_dist_xx/(2*exp(2*ll)) - sq_dist_tt/(2*exp(2*llt))).*sq_dist_xx;
elseif index == 2	% dK/dlt
	sq_dist_xx = square_dist(x, x);
	sq_dist_tt = square_dist(t, t);
	C = exp(2*lsf - 2*llt - sq_dist_xx/(2*exp(2*ll)) - sq_dist_tt/(2*exp(2*llt))).*sq_dist_tt;
elseif index == 3	% dK/dsf
	sq_dist_xx = square_dist(x, x);
	sq_dist_tt = square_dist(t, t);
	C = 2*exp(lsf - sq_dist_xx/(2*exp(2*ll)) - sq_dist_tt/(2*exp(2*llt)));
else				% dK/dsn
	C = 2*exp(2*loghyper(4))*eye(size(x, 1));
end