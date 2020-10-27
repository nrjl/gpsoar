function C = dsquare_exp_dx(a, b, k, loghyper)
% loghyper: log([l, sigma_n, sigma_f, sigma_n])
n = size(b, 1);		% Number of observations
m = size(a, 1);		% Number of test points

C = -(repmat(a(:,k),[1,n]) - repmat(b(:,k)',[m,1]))/exp(2*loghyper(1)).*...
	exp(2*loghyper(2)-0.5/exp(2*loghyper(1))*square_dist(a, b));