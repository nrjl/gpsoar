function C = square_exp_noniso(a, b, loghyper, varargin)

if nargin == 5
	b = loghyper;
	loghyper = varargin{2};
end
% loghyper: log([l1, L2, ... ln, sigma_n, sigma_f, sigma_n])
l = loghyper(1:end-3);
C = exp(2*loghyper(2)-0.5/exp(2*loghyper(1))*square_dist(a, b, l));