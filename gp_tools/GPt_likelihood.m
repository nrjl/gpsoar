function [nlml, dnlml] = GPt_likelihood(x, t, y, cov_funs, loghyper)
%--------------------------------------------------------------------------
%
% FUNCTION:		GPt_likelihood
%
% PURPOSE:		Calculate Gaussian Process log marginal likelihood and its
%				derivative (TEMPORAL VERSION)
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
%				loghypert- log of temporal hyperparameters [n_hypert×1]
%
% OUTPUTS:		nlml	- negative log marginal likelihood [1×1]
%				dnlml	- derivative of nlml wrt hyperparams 
%							[(n_hyper+n_hypert)×1]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2009
%
% MODIFIED:     November 2009
%
% See also:		GP_likelihood, GPt_likelihoodn, GPt_predict
%--------------------------------------------------------------------------

% loghyper = log([length_scale, lag_scale, sigma_f, sigma_n]);
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
Kxt = cov_funs{1}(x, t, x, t, loghyper);

% loghyper(end) = max([loghyper(end), log(1)]);

K_mat = Kxt + exp(2*loghyper(end))*eye(length(x));

% Cholesky decomposition:
% K*A = y
% K = L*L'
% Use Matlab Cholesky decomp to get L (lower trangular)
% Lz = y	=>	z = L\y
% L'x = z	=>	A = L'\z
%			=>	A = L'\(L\y)
adanoise=1e-6;
while(true)
	try
		L = chol(K_mat, 'lower');
		break;
	catch
		fprintf('!');
		K_mat = K_mat + adanoise*eye(size(K_mat));
		adanoise = adanoise*10;
	end
end

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

nlml = 0.5*y'*alpha + sum(log(diag(L))) + 0.5*n_x*log(2*pi)+hyper_prior;

if nargout == 2
	n_hyper = numel(loghyper);
	dnlml = zeros(n_hyper, 1);
	inner = alpha*alpha' - K_inv;
	
	if (numel(cov_funs) > 1) && ~isempty(cov_funs{2})
		for i = 1:n_hyper
			% Normally use trace(inner*dcov) but for symmetrical, this is
			% equivalent to sum(sum(inner.*dcov))
			dnlml(i) = -0.5*sum(sum(inner.*cov_funs{2}(x, t, i, loghyper))); % some comment, you dont need to know
		end
	else
		for i = 1:n_hyper
			newhyper = loghyper;
			newhyper(i) = newhyper(i) + 1e-6;
			nlml2 = GPt_likelihood(x, t, y, cov_funs, [newhyper; lh0; sigma_hyper]);
			dnlml(i) = (nlml2-nlml)/(1e-6);
		end
	end
end