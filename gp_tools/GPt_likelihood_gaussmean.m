function [nlml] = GPt_likelihood_gaussmean(x, t, y, cov_funs, mean_fun, loghyper, n_hyper)
%--------------------------------------------------------------------------
%
% FUNCTION:		GPt_likelihood_gaussmean
%
% PURPOSE:		Calculate Gaussian Process log marginal likelihood and its
%				derivative (TEMPORAL VERSION)
%               
% SYNTAX:		[nlml, dnlml] = ...
%					 GPt_likelihood_gaussmean(x, t, y, cov_funs, mean_fun,
%					 loghyper, n_hyper)
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
% CREATED:		November 2009
%
% MODIFIED:     November 2009
%
% See also:		GP_likelihood, GPt_likelihoodn, GPt_predict
%--------------------------------------------------------------------------

% sq_expt: loghyper = log([length_scale, lag_scale, sigma_f, (sigma_n)]);
% sq_exp : loghyper = log([length_scale, sigma_f, (sigma_n)]);
n_x = size(x, 1);
n_nonmeanhyper = sum(n_hyper(1:end-1));

Kxt = cov_funs{1}(x, t, x, t, loghyper(1:n_nonmeanhyper));

H = mean_fun(x);
order_h = size(H, 1);

b = (loghyper((n_nonmeanhyper+1):(n_nonmeanhyper+order_h)))';
% B = diag(abs((loghyper((n_nonmeanhyper+1+order_h):(n_nonmeanhyper+2*order_h)))), 0);
B = diag((loghyper((n_nonmeanhyper+1+order_h):(n_nonmeanhyper+2*order_h))).^2, 0);
	
K_mat = Kxt + exp(2*loghyper(n_hyper(1)))*eye(length(x));
K_mod = K_mat + H'*B*H;

% Cholesky decomposition:
% K*A = y
% K = L*L'
% Use Matlab Cholesky decomp to get L (lower trangular)
% Lz = y	=>	z = L\y
% L'x = z	=>	A = L'\z
%			=>	A = L'\(L\y)

L_mod = chol(K_mod, 'lower');
K_mod_inv = (L_mod'\(L_mod\eye(n_x)));

% alpha = K_inv*y;

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

% --- Supplied priors on beta ~ N(b, B) --- %
nlml = 0.5*(H'*b - y)'*K_mod_inv*(H'*b - y) + sum(log(diag(L_mod))) + ...
	0.5*n_x*log(2*pi);

%% SCRAP

% ----- DNLML OUTPUT ----- %
% n_hyper = numel(loghyper);
% 
% dnlml = zeros(n_hyper, 1);
% inner = alpha*alpha' - K_inv;
% 
% for i = 1:n_hyper
% 	dnlml(i) = 0.5*trace(inner*cov_funs{2}(x, t, i, loghyper));
% end


% --- Normal version --- %
% nlml = 0.5*y'*alpha + sum(log(diag(L))) + 0.5*n_x*log(2*pi);

% --- Vague prior (inv(B) ~ [0]) --- %
% A = H*K_inv*H';
% C = K_inv*H'*(A\H)*K_inv;
% nlml = 0.5*y'*alpha - 0.5*y'*C*y + sum(log(diag(L))) +0.5*log(det(A)) +
% (n_x-m)/2*log(2*pi);