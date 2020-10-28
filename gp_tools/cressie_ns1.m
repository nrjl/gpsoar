function K = cressie_ns1(x1, t1, x2, t2, hyper, varargin)
% hyper = [sigma, lt, l]
if nargin == 5
	h2 = square_dist(x1, x2);
	u2 = square_dist(t1, t2);
else
	h2 = varargin{1};
	u2 = varargin{2};
end

a2_u2 = (.5*hyper(2)^(-2)*u2+1);
d = size(x1, 2);

K = hyper(1).^2./(a2_u2.^(d/2)).*exp(-0.5*hyper(3).^(-2)*h2./a2_u2);