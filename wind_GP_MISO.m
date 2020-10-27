function [f_star, V, log_marginal] = wind_GP_MISO(x, y, x_star, cov_fun, loghyper)
%% GP Procedure
%
% x = input training points
% y = output training points
%
% x_star = test point input
% f_star = test point output
%
% loghyper = log(hyperparameters);
%		   = log([length_scale, sigma_f, sigma_n]) for SE
% cov_fun = handle to covariance function
%

if nargin == 3
	length_scale = 1;
	sigma_f = 0.5;
	sigma_n = 0.5;
	cov_fun = @square_exp;
	loghyper = log([length_scale, sigma_f, sigma_n]);
end


n_obs = numel(y);

K = cov_fun(x, x, loghyper);				% Training points covariances
K_star = cov_fun(x_star, x_star, loghyper);	% Test point covariances
k_star = cov_fun(x, x_star, loghyper);		% Related covariances


K_mat = K + exp(2*loghyper(3))*eye(length(x));

% Cholesky decomposition:
% K*A = y
% K = L*L'
% Use Matlab Cholesky decomp to get L (lower trangular)
% Lz = y	=>	z = L\y
% L'x = z	=>	A = L'\z
%			=>	A = L'\(L\y)

L = chol(K_mat, 'lower');
K_inv_y = (L'\(L\y));
f_star = k_star'*K_inv_y;


if nargout == 2
	V = K_star - k_star'*(L'\(L\k_star));          % Full covariance
elseif nargout == 3
	V = K_star - k_star'*(L'\(L\k_star));          % Full covariance
	log_marginal = -0.5*(y'*K_inv_y - log(det(K + exp(2*loghyper(3))*...
		eye(length(x)))) - n_obs*log(2*pi));
end


%% SCRAP
% L = chol(K_mat);    Ri = R^-1;
% K_inv = Ri*Ri';

% K_inv = inv(K_mat);						% Inversion
% f_star = k_star'*K_inv*y;					% Outputs
