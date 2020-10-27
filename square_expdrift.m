function Kxd = square_expdrift(x1, t1, x2, t2, loghyper)
% loghyper = log([l, sigma_f, wx, wy, wz, sigma_n]);

x_dim = size(x1, 2);

n_l = size(loghyper, 2) - x_dim - 2;
W = loghyper((n_l+2):(n_l+1+x_dim));

dd = square_distdrift(x1, t1, x2, t2, W);
Kxd = exp(2*loghyper(2)-0.5/exp(2*loghyper(1))*dd);
