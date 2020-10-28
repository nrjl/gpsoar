function K = d_cressie_ns1(x, t, index, hyper)
% hyper = [sigma, lt, l (, log(sn))]

h2 = square_dist(x, x);
u2 = square_dist(t, t);

a2_u2 = (hyper(2)^(-2)*u2+1);
d = size(x, 2);

K_temp = exp(-hyper(3).^(-2)*h2./a2_u2);


switch index
	case 1
		K = 2*hyper(1).*(a2_u2.^(-d/2)).*K_temp;
	case 2
% 		K = -(u2*d*hyper(2)^(-3) + hyper(2).^2*h2*hyper(3).^(-2))*a2_u2.^(-3*d/2).*K_temp;
		K = hyper(1)^2*hyper(2)^(-3)*u2.*(d + 2*hyper(3)^(-2)*h2.*((a2_u2).^(-1))).*((a2_u2).^(-1-d/2)).*K_temp;
	case 3
		K = 2*hyper(1).^2*h2*hyper(3)^(-3).*a2_u2.^(-3*d/2).*K_temp;
	case 4
		K = 2*exp(2*hyper(4))*eye(size(x, 1));
end