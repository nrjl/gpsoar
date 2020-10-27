function K = cressie_ns5(x1, t1, x2, t2, hyper, varargin)
% hyper = [sf, lt, l, c]

if nargin == 5
	h2 = square_dist(x1, x2);
	u2 = square_dist(t1, t2);
else
	h2 = varargin{1};
	u2 = varargin{2};
end

lti2 = .5*hyper(2)^(-2);
li2 = .5*hyper(3)^(-2);

K = hyper(1)^2*exp(-lti2*u2 - li2*h2 - abs(hyper(4))*(u2*lti2).*(h2*li2));