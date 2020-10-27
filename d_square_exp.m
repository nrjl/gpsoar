function C = d_square_exp(x, index, loghyper)
% Different from Rasmussen!!
% loghyper: log([l, sigma_f, sigma_n])
ll = loghyper(1);
lsf = loghyper(2);

if index == 1		% dK/dl
	sq_dist_xx = square_dist(x, x);
	C = exp(2*lsf - 2*ll - sq_dist_xx/(2*exp(2*ll))).*sq_dist_xx;
elseif index == 2	% dK/dsf
	C = 2*exp(2*lsf - square_dist(x, x)/(2*exp(2*ll)));
else				% dK/dsn
	C = 2*exp(2*loghyper(3))*eye(size(x, 1));
end