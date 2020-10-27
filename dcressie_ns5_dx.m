function K = dcressie_ns5_dx(x1, t1, x2, t2, k, hyper)
% hyper = [sf, lt, l, c]
n = size(x2, 1);		% Number of observations
m = size(x1, 1);		% Number of test points

h2 = square_dist(x1, x2);
u2 = square_dist(t1, t2);

lti2 = hyper(2)^(-2);
li2 = hyper(3)^(-2);


K = -(repmat(x1(:,k),[1,n]) - repmat(x2(:,k)',[m,1])).*(li2+hyper(4)*u2).*...
	hyper(1)^2*exp(-lti2*u2 - li2*h2 - abs(hyper(4))*li2*lti2*u2.*h2);
