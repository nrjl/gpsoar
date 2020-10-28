function [nlml] = GPtp_likelihood(x, t, y, cov_funs, loghyper)
%--------------------------------------------------------------------------
%
% FUNCTION:		GPtp_likelihood
%
% PURPOSE:		Calculate Gaussian Process log marginal likelihood and its
%				derivative (TEMPORAL & PROJECTION VERSION)
%               
% SYNTAX:		[nlml, dnlml] = ...
%					GPt_likelihood(x, t, y, cov_funs, loghyper, loghypert)
%
% INPUTS:		x		- input training points [n×dim]
%				t		- input times [n×1]
%				y		- output training points [n×1]
%				cov_funs- handles to covariance function and its derivative
%							in structure form:	{cov_fun, tcov_fun, dcov_fun} 
%				loghyper - log of spatial hyperparameters [n_hyper×1]
%
% OUTPUTS:		nlml	- negative log marginal likelihood [1×1]
%				dnlml	- derivative of nlml wrt hyperparams 
%							[(n_hyper+n_hypert)×1]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		March 2010
%
% MODIFIED:     March 2010
%
% See also:		GP_likelihood, GPt_likelihood, GPtp_predict
%--------------------------------------------------------------------------

% loghyper = log([lx, lt, lp, sigma_f, sigma_n]);
n_x = size(x, 1);
t_test = t(round(numel(t)/2));

d_xx = square_dist(x, x);
d_tt = square_dist(t, t);

x_prime = x + diag(t_test-t,0)*y;
d_pp = square_dist(x_prime, x_prime);

Kxp = projected_square_expt(x, t, y, x, t, loghyper(1:end-1), d_xx, d_tt, d_pp);

K_mat = Kxp + exp(2*loghyper(end))*eye(length(x));

% Cholesky decomposition:
% K*A = y
% K = L*L'
% Use Matlab Cholesky decomp to get L (lower trangular)
% Lz = y	=>	z = L\y
% L'x = z	=>	A = L'\z
%			=>	A = L'\(L\y)

L = chol(K_mat, 'lower');
K_inv = (L'\(L\eye(n_x)));

alpha = K_inv*y;

% Maths for det(K)
% Standard equation is:
% nlml = 0.5*y'*alpha + 0.5*log(det(K_mat)) + 0.5*n_x*log(2*pi);
% However, middle term is badly conditioned, resulting in bad results.
% Then, use some tricky maths:
% 0.5*log|K| = 0.5*log(|L| * |L'|)
%			 = log(|L|)
%			 = log(exp(tr(log(L))))
%			 = tr(log(L))
%			 = sum(log(diag(L))

nlml = mean(diag(0.5*y'*alpha)) + sum(log(diag(L))) + 0.5*n_x*log(2*pi);

% n_hyper = numel(loghyper);
% 
% dnlml = zeros(n_hyper, 1);
% inner = alpha*alpha' - K_inv;
% 
% for i = 1:n_hyper
% 	% Normally use trace(inner*dcov) but for symmetrical, this is 
% 	% equivalent to sum(sum(inner.*dcov))
% 	dnlml(i) = -0.5*sum(sum(inner.*cov_funs{2}(x, t, i, loghyper))); % some comment, you dont need to know
% end