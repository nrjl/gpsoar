function Kxd = square_exptdrift(x1, t1, x2, t2, loghyper)
% loghyper = (log([l, lt, ld, sigma_f]), K, wx, wy, wz, log(sigma_n)]);

x_dim = size(x1, 2);

n_l = 3;
W = loghyper((n_l+3):(n_l+2+x_dim));
lxi2 = exp(-2*loghyper(1));
lti2 = exp(-2*loghyper(2));
ldi2 = exp(-2*loghyper(3));

dd = square_distdrift(x1, t1, x2, t2, W);
xx = square_dist(x1,x2);
tt = square_dist(t1,t2);

KK = 1/(1+exp(-loghyper(5)));

% Kxd = exp(2*loghyper(4)-0.5*lti2*tt).*(KK*exp(-0.5*lxi2*xx) + (1-KK)*exp(-0.5*ldi2*dd));
Kxd = exp(2*loghyper(4)-0.5*lti2*tt -0.5*KK*lxi2*xx - 0.5*(1-KK)*ldi2*dd);
