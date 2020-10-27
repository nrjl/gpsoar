function C = square_exp(a, b, loghyper, varargin)
% loghyper: log([l1, l2, ... ln, sigma_f, sigma_n])
% can supply square distance, always last argument (ie 4 or 6)


if rem(nargin, 2)		%nargin = 3 or 5 -> calculate sqdist
	if nargin == 5
		b = loghyper;
		loghyper = varargin{2};
	end
	
	if numel(loghyper) > 3
		d = size(a,2);	% number of dimensions
		l = exp(loghyper(1:d));
		C = exp(2*loghyper(d+1)-0.5*square_dist(a, b, l));
	else
		C = exp(2*loghyper(2)-0.5/exp(2*loghyper(1))*square_dist(a, b));
	end
else					% nargin = 4, 6 -> sq_dist supplied
	sqdist = varargin{(nargin-3)};
	if numel(loghyper) > 3
		ll = 0.5; d = size(a,2);
	else
		ll = 0.5/exp(2*loghyper(1)); d = 1;
	end

	C = exp(2*loghyper(d+1)- ll*sqdist);
end

% Old version:
% 
% if numel(loghyper) > 3
% 	d = size(a,2);	% number of dimensions
% 	l = exp(loghyper(1:d));
% 	C = exp(2*loghyper(d+1)-0.5*square_dist(a, b, l));
% else
% 	C = exp(2*loghyper(2)-0.5/exp(2*loghyper(1))*square_dist(a, b));
% end
