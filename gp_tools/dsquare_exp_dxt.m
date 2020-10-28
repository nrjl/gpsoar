function C = dsquare_exp_dxt(x1, t1, x2, t2, k, loghyper)
% loghyper: log([l, sigma_n, sigma_f, sigma_n])
n = size(x2, 1);		% Number of observations
m = size(x1, 1);		% Number of test points

C = -(repmat(x1(:,k),[1,n]) - repmat(x2(:,k)',[m,1]))/exp(2*loghyper(1)).*...
	exp(2*loghyper(3)-0.5/exp(2*loghyper(1))*square_dist(x1, x2) ...
					 -0.5/exp(2*loghyper(2))*square_dist(t1, t2));