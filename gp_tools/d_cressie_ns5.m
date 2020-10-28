function K = d_cressie_ns5(x, t, index, hyper)
% hyper = [sf, lt, l, c, log(sn)]

h2 = square_dist(x, x);
u2 = square_dist(t, t);

lti2 = .5*hyper(2)^(-2);
li2 = .5*hyper(3)^(-2);

K_temp = exp(-lti2*u2 - li2*h2 - abs(hyper(4))*(u2*lti2).*(h2*li2));

switch index
	case 1
		K = 2*hyper(1)*K_temp;
	case 2
		K = hyper(1).^2*(2*hyper(2)^(-3))*(u2 + abs(hyper(4))*li2*u2.*h2).*K_temp;
	case 3
		K = hyper(1).^2*(2*hyper(3)^(-3))*(h2 + abs(hyper(4))*lti2*u2.*h2).*K_temp;
	case 4
		K = hyper(1).^2*-(u2*lti2).*(h2*li2).*K_temp;
	case 5
		K = 2*exp(2*hyper(5))*eye(size(x, 1));
end