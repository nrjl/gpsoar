function [f_star, V] = GPtp_predict(X, t, Y, X_star, t_star, loghyper, K_inv)
%--------------------------------------------------------------------------
%
% FUNCTION:		GPtp_predict
%
% PURPOSE:		Make a Gaussian Process prediction and return mean and
%				covariance at specified target points
%               
% SYNTAX:		[f_star, V] =  GPtp_predict(X, t, Y, X_star, t_star,
%					loghyper, K_inv)
%
% INPUTS:		X		- input training points
%				Y		- output training points
%				X_star	- test point input
%				loghyper- log of hyperparameters
%
% OUTPUTS:		f_star	- test point output
%				V		- full covariance
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		March 2010
%
% MODIFIED:     March 2010
%
% See also:		
%--------------------------------------------------------------------------

% Get inv(K)
if nargin ~= 5
	
	n_obs = size(X, 1);
	
	d_xx = square_dist(X, X);
	d_tt = square_dist(t, t);
	
	X_prime = X + diag(t_star-t,0)*Y;
	d_pp = square_dist(X_prime, X_prime);
	
	K = projected_square_expt(X, t, Y, X, t, loghyper, d_xx, d_tt, d_pp);
	
	K_mat = K + exp(2*loghyper(end))*eye(n_obs); % Add noise
	
	% Cholesky decomposition:
	% K*A = Y
	% K = L*L'
	% Use Matlab Cholesky decomp to get L (lower trangular)
	% Lz = Y	=>	z = L\Y
	% L'X = z	=>	A = L'\z
	%			=>	A = L'\(L\Y)
	
	L = chol(K_mat, 'lower');
	K_inv = (L'\(L\eye(n_obs)));	
end


% OUTPUT
if isempty(X_star)
	f_star = K_inv;
else
	d_xxs = square_dist(X, X_star);
	d_tts = square_dist(t, t_star*ones(size(X_star, 1),1));
	
	d_pps = square_dist(X_prime, X_star);
	
	k_star = projected_square_expt(X, t, Y, X_star, t_star, loghyper, d_xxs, d_tts, d_pps)';
	
	f_star = k_star*K_inv*Y;

	if nargout == 2
		% Limited (self-) covariances of test points
		n_test = size(X_star, 1);
		V = exp(2*loghyper(end))*ones([n_test, 1]) - ...
			reshape(diag(k_star*K_inv*k_star'), [n_test, 1]);
	end
end

% K_star = cov_fun(X_star, X_star, loghyper);	% Test point covariances
% V = K_star - k_star'*K_inv*k_star;		% Full covariance