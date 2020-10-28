function [nlml, dnlml] = GP_likelihood(x, y, cov_funs, loghyper)
%--------------------------------------------------------------------------
%
% FUNCTION:		GP_likelihood
%
% PURPOSE:		Calculate Gaussian Process log marginal likelihood and its
%				derivative
%               
% SYNTAX:		[nlml, dnlml] = GP_likelihood(x, y, cov_funs, loghyper)
%
% INPUTS:		x		- input training points [n×dim]
%				y		- output training points [n×1]
%				cov_funs- handles to covariance function and its derivative
%							in structure form:	{cov_fun, dcov_fun} 
%				loghyper- log of hyperparameters [n_hyper×1]
%
% OUTPUTS:		nlml	- negative log marginal likelihood [1×1]
%				dnlml	- derivative of nlml wrt hyperparams [n_hyper×1]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		June 2009
%
% MODIFIED:     July 2009
%
% See also:		GP_predict
%--------------------------------------------------------------------------

if size(loghyper, 1) > 1
	lh0 = loghyper(2,:);
	sigma_hyper = loghyper(3,:);
	loghyper = loghyper(1,:);
	dtheta = exp(lh0) - exp(loghyper);
	hyper_prior = dtheta.*exp(-2*sigma_hyper)*dtheta';
else
	hyper_prior = 0;
end

n_x = size(x, 1);
K = cov_funs{1}(x, x, loghyper);			% Training point covariances

K_mat = K + exp(2*loghyper(end))*eye(length(x));

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

nlml = 0.5*y'*alpha + sum(log(diag(L))) + 0.5*n_x*log(2*pi)+ hyper_prior;

if nargout == 2
	n_hyper = numel(loghyper);
	dnlml = zeros(n_hyper, 1);
	inner = alpha*alpha' - K_inv;
	
	if numel(cov_funs) > 1
		for i = 1:n_hyper
			dnlml(i) = -0.5*trace(inner*cov_funs{2}(x, i, loghyper));
		end
	else
		for i = 1:n_hyper
			newhyper = loghyper;
			newhyper(i) = newhyper(i) + 1e-8;
			nlml2 = GP_likelihood(x, y, cov_funs, newhyper);
			dnlml(i) = (nlml2-nlml)/1e-8;
		end
	end
end