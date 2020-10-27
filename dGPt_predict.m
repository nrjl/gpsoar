function [df_star] = dGPt_predict(x, t, y, x_star, t_star, k, cov_fun, loghyper, K_inv)
%--------------------------------------------------------------------------
%
% FUNCTION:		dGPt_predict
%
% PURPOSE:		Make a Gaussian Process prediction to estimate mean
%				derivative with respect to kth input dimension at specified
%				points (x_star, t_star) from observations (x, t, y).
%               
% SYNTAX:		[df_star] = 
%				 dGP_predict(x, t, y, x_star, t_star, k, cov_fun, loghyper)
%
%				To just provide inv(K)*y product use empty test point set:
%				[Ki] = dGPt_predict(x, t, y, [], [], cov_fun, loghyper)
%
%				To solve for a known Ki = inv(K):
%				[df] = dGP_predict(x, t, y, x_star, t_star, k, cov_fun, ...
%					loghyper, Ki)
%
% INPUTS:		x		- input training points
%				t		- input training times
%				y		- output training points
%				x_star	- test point input
%				t_star	- test point times
%				k		- input dimension partial derivative index(es)
%				cov_fun - handle to covariance function(s). If Ki is 
%	specified, then the user only needs to supply the derivative covariance 
%	function. If Ki is not supplied, then the original and derivative must 
%	be supplied ({@cov, @dcov})
%				loghyper- log of hyperparameters
%
% OUTPUTS:		f_star	- test point output
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		September 2009
%
% MODIFIED:     July 2010
%
% See also:		GPt_predict, GPt_likelikood, dGP_predict
%--------------------------------------------------------------------------

% Get inv(K)
if nargin ~= 9
	n_obs = size(y, 1);
	
	K = cov_fun{1}(x, t, x, t, loghyper);			% Training points covariances
	K_mat = K + exp(2*loghyper(3))*eye(length(x)); % Add noise
	
	% Cholesky decomposition:
	% K*A = y
	% K = L*L'
	% Use Matlab Cholesky decomp to get L (lower triangular)
	% Lz = y	=>	z = L\y
	% L'x = z	=>	A = L'\z
	%			=>	A = L'\(L\y)
	
	L = chol(K_mat, 'lower');
	K_inv = (L'\(L\eye(n_obs)));
	cov_fun = cov_fun{2};
end


% OUTPUT
if isempty(x_star)
	df_star = K_inv;
else
	if numel(k) > 1
		df_star = zeros(numel(k), size(x, 2), size(x_star, 1));
		for i = 1:numel(k)
			k_star = cov_fun(x_star, t_star, x, t, k(i), loghyper);	% Related covariances
			df_star(i,:,:) = (k_star*K_inv*y)';
		end
	else
		k_star = cov_fun(x_star, t_star, x, t, k, loghyper);	% Related covariances
		df_star = k_star*K_inv*y;
	end
end