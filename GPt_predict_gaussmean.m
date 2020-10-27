function [f_star, V] = GPt_predict_gaussmean(x, t, y, x_star, t_star, cov_fun, mean_fun, loghyper, n_hyper, K_inv)
%--------------------------------------------------------------------------
%
% FUNCTION:		GPt_predict_gaussmean
%
% PURPOSE:		Make a Gaussian Process prediction with a specified mean 
%				function with known prior and return mean and covariance
%				at specified target points
%               
% SYNTAX:		[f_star] = GPt_predict_gaussmean(x, t, y, x_star, t_star, ...
%					cov_fun, mean_fun, loghyper, n_hyper)
%				[f_star, V] = GPt_predict_gaussmean(x, t, y, x_star, t_star, ...
%					cov_fun, mean_fun, loghyper, n_hyper)
%
%				To just provide inv(K) product use empty test point set:
%				[Ki] = GPt_predict_gaussmean(x, t, y, [], [], cov_fun, loghyper, n_hyper)
%
%				To solve for a known Ki = inv(K):
%				[f_star, V] = GPt_predict_gaussmean(x, t, y, x_star, t_star, ...
%					cov_fun, mean_fun, loghyper, n_hyper, Ki)
%
% INPUTS:		x		- input training points
%				t		- input training time vector
%				y		- output training points
%				x_star	- test points
%				t_star	- test times
%				cov_fun - handle to covariance function
%				mean_fun- handle to (explicit) mean function
%				loghyper- log of hyperparameters
%				n_hyper - number of hyperparameters per function
%
% OUTPUTS:		f_star	- test point output
%				V		- full covariance
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2009
%
% MODIFIED:     November 2009
%
% See also:		GPt_predict
%--------------------------------------------------------------------------

n_nonmeanhyper = sum(n_hyper(1:end-1));

% Get inv(K)
if nargin ~= 10
	
	n_obs = size(x, 1);
	
	K = cov_fun(x, t, x, t, loghyper(1:n_nonmeanhyper));			% Training points covariances
	K_mat = K + exp(2*loghyper(n_hyper(1)))*eye(size(K)); % Add noise
	
	% Cholesky decomposition:
	% K*A = y
	% K = L*L'
	% Use Matlab Cholesky decomp to get L (lower trangular)
	% Lz = y	=>	z = L\y
	% L'x = z	=>	A = L'\z
	%			=>	A = L'\(L\y)
	
	L = chol(K_mat, 'lower');
	K_inv = (L'\(L\eye(n_obs)));	
end


% OUTPUT
if isempty(x_star)
	f_star = K_inv;
else
	n_nonmeanhyper = sum(n_hyper(1:end-1));
	k_star = cov_fun(x_star, t_star, x, t, loghyper(1:n_nonmeanhyper));		% Related covariances
	
	H = mean_fun(x);
	H_star = mean_fun(x_star);
	order_h = size(H, 1);
	
	b = (loghyper((n_nonmeanhyper+1):(n_nonmeanhyper+order_h)))';
% 	B_inv = diag(1./abs(loghyper((n_nonmeanhyper+1+order_h):(n_nonmeanhyper+2*order_h))), 0);
	B_inv = diag(1./(loghyper((n_nonmeanhyper+1+order_h):(n_nonmeanhyper+2*order_h))).^2, 0);
% 	B_inv = diag(1./(loghyper((5+order_h):(4+2*order_h))), 0);
	
	Rt = H_star' - k_star*K_inv*H';
	beta_bar = (B_inv + H*K_inv*H')\(H*K_inv*y + B_inv*b);
	
	f_star = k_star*K_inv*y + Rt*beta_bar;
	% f_star = H_star'*beta_bar + k_star*K_inv*(y - H'*beta_bar);


	if nargout == 2
		% Limited (self-) covariances of test points
		n_test = size(x_star, 1);
		% V = exp(2*loghyper(3))*ones([n_test, 1]) - ...
		% 	reshape(diag(k_star*K_inv*k_star'), [n_test, 1]) + ...
		% 	reshape(diag(R'*((H*K_inv*H')\R)), [n_test, 1]);
		
% 		V = exp(2*loghyper(3))*ones([n_test, 1]) - ...
% 			reshape(diag(k_star*K_inv*k_star'), [n_test, 1]) + ...
% 			reshape(diag(Rt*((B_inv + H*K_inv*H')\Rt')), [n_test, 1]);	

		% This is the slow (full) version which doesn't use the sq_exp hack
		% to quickly calculate the self covariances of the test points. 
		V = reshape(diag(cov_fun(x_star, t_star, x_star, t_star, loghyper)) - ...
			diag(k_star*K_inv*k_star') + diag(Rt*((B_inv + H*K_inv*H')\Rt')), [n_test, 1]);
			
	end
end

% K_star = cov_fun(x_star, x_star, loghyper);	% Test point covariances
% V = K_star - k_star'*K_inv*k_star;		% Full covariance