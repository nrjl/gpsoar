function [f_star, V] = GPt_predictn_gaussmean(x, t, y, x_star, t_star, cov_fun, mean_fun, loghyper, n_hyper, K_inv)
%--------------------------------------------------------------------------
%
% FUNCTION:		GPt_predictn_gaussmean
%
% PURPOSE:		Make a Gaussian Process prediction with a specified mean 
%				function with unknown prior and return mean and covariance
%				at specified target points
%               
% SYNTAX:		[nlml, dnlml] = ...
%					GP_likelihoodn(x, t, y1, y2, ... yn, cov_funs, mean_fun, loghyper, n_hyper)
%
% INPUTS:		x		- input training points [i×dim]
%				yk		- output training points k [n×1]
%				cov_funs- handles to covariance function and its derivative
%							in structure form:	{cov_fun, dcov_fun}
%				mean_fun- handle to mean basis functions
%				loghyper- log of hyperparameters [n_hyper×1]
%				n_hyper - number of hyperparameters per function
%
% OUTPUTS:		nlml	- negative log marginal likelihood [1×1]
%				dnlml	- derivative of nlml wrt hyperparams [n_hyper×1]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2009
%
% MODIFIED:     August 2009
%
% See also:		GP_likelihood, GP_predict
%--------------------------------------------------------------------------

n_y = size(y, 2);
order_h = round(n_hyper(end)/(2*n_y));
n_nonmeanhyper = sum(n_hyper(1:end-1));

f_star = zeros(size(x_star, 1), n_y);
V = f_star(:,1);

if nargin < 10
	[K_inv] = GPt_predict_gaussmean(x, t, y(:,1), [], [], cov_fun, mean_fun, loghyper, n_hyper);
end

for i = 1:n_y
	loghyper_temp = loghyper([1:n_nonmeanhyper, ((n_nonmeanhyper+1+(i-1)*2*order_h):(n_nonmeanhyper+i*order_h*2))]);
	[f_star(:,i), V(:,1)] = GPt_predict_gaussmean(x, t, y(:,i), x_star, t_star,...
		cov_fun, mean_fun, loghyper_temp, n_hyper, K_inv);
end