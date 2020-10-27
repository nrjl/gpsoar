function C = dsquare_exptdrift_dxt(x1, t1, x2, t2, k, loghyper)
% loghyper = (log([l, lt, ld, sigma_f]), K, wx, wy, wz, log(sigma_n)]);
n = size(x2, 1);		% Number of observations
m = size(x1, 1);		% Number of test points

x_dim = size(x1, 2);
n_l = 3;
W = loghyper((n_l+3):(n_l+2+x_dim));
KK = 1/(1+exp(-loghyper(5)));

ddx = repmat(x1(:,k),[1,n]) - repmat(x2(:,k)',[m,1]);

C = -( KK*ddx*exp(-2*loghyper(1)) + (1-KK)*exp(-2*loghyper(3))*(ddx + ...
	W(k)*(repmat(t1,[1,n]) - repmat(t2',[m,1])))).*...
	square_exptdrift(x1, t1, x2, t2, loghyper);