function [f_star, V] = GPt_predict(x, t, y, x_star, t_star, cov_fun, loghyper, K_inv)
%--------------------------------------------------------------------------
%
% FUNCTION:		GPt_predict
%
% PURPOSE:		Make a Gaussian Process prediction and return mean and
%				covariance at specified target points
%               
% SYNTAX:		[f_star] = GP_predict(x, y, x_star, cov_fun, loghyper)
%				[f_star, V] = GP_predict(x, y, x_star, cov_fun, loghyper)
%
%				To just provide inv(K)*y product use empty test point set:
%				[Ki] = GP_predict(x, y, [], cov_fun, loghyper)
%
%				To solve for a known Ki = inv(K):
%				[f_star] = GP_predict(x, y, x_star, cov_fun, loghyper, Ki)
%				[f_star, V] = GP_predict(x,y,x_star, cov_fun, loghyper, Ki)
%
% INPUTS:		x		- input training points
%				y		- output training points
%				x_star	- test point input
%				cov_fun - handle to covariance function
%				loghyper- log of hyperparameters
%
% OUTPUTS:		f_star	- test point output
%				V		- full covariance
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		June 2009
%
% MODIFIED:     July 2010
%
% See also:		
%--------------------------------------------------------------------------

% Get inv(K)
if nargin ~= 8
	
	n_obs = size(x, 1);
	
	K = cov_fun(x, t, x, t, loghyper);			% Training points covariances
	K_mat = K + exp(2*loghyper(end))*eye(n_obs); % Add noise
	
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
	K_inv = (L'\(L\eye(n_obs)));	
end


% OUTPUT
if isempty(x_star)
	f_star = K_inv;
else
	k_star = cov_fun(x_star, t_star, x, t, loghyper);		% Related covariances
	f_star = k_star*K_inv*y;

	if nargout == 2
		% Limited (self-) covariances of test points
		n_test = size(x_star, 1);
		V = exp(2*loghyper(4))*ones([n_test, 1]) - ...
			reshape(diag(k_star*K_inv*k_star'), [n_test, 1]);
	end
end

% K_star = cov_fun(x_star, x_star, loghyper);	% Test point covariances
% V = K_star - k_star'*K_inv*k_star;		% Full covariance