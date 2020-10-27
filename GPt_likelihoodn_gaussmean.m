function [nlml] = GPt_likelihoodn_gaussmean(x, t, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		GPt_likelihoodn_gaussmean
%
% PURPOSE:		Calculate GP negative log marginal likelihood and 
%				derivative for n inputs with gaussian mean function
%               
% SYNTAX:		[nlml, dnlml] = GPt_likelihoodn_gaussmean(x, t, y1, y2, ... yn, ...
%					cov_funs, mean_fun, loghyper, n_hyper)
%
% INPUTS:		x		- input training points [i×dim]
%				yk		- output training points k [n×1]
%				cov_funs- handles to covariance function and its derivative
%							in structure form:	{cov_fun, dcov_fun}
%				mean_fun- handle to mean basis functions
%				loghyper- log of hyperparameters [n_hyper×1]
%				n_hyper - Number of hyperparameters per function. Note:
% n_hyper is a list of the number of hyperparameters of the corresponding
% covariance function. Any arguments for the mean function must be appended
% AFTER the other covariance functions, and should still have a 
% corresponding entry in the n_hyper vector, but will not have a 
% fun_handles entry (i.e numel(fun_handles) = numel(n_hyper) or 
% numel(fun_handles) = numel(n_hyper)-1.
%
% Example:  For a sum of two covariance functions with 4 HPs each, and a 
% mean function with 6 HPs the corresponding n_hyper is [4 4 6].
%
% OUTPUTS:		nlml	- negative log marginal likelihood [1×1]
%				dnlml	- derivative of nlml wrt hyperparams [n_hyper×1]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2009
%
% MODIFIED:     December 2009
%
% See also:		GP_likelihood, GPt_likelihood, GPt_likelihood_gaussmean
%--------------------------------------------------------------------------

if nargin < 7
	error('GPNick:GP_likelihoodn:n_inputs', 'Too few input arguments');
else
	cov_funs = varargin{end-3};
	mean_fun = varargin{end-2};
	loghyper = varargin{end-1};
	n_hyper = varargin{end};
	n_y = nargin - 6;
end

nlml = 0;
% dnlml = 0;
order_h = round(n_hyper(end)/(2*n_y));
n_nonmeanhyper = sum(n_hyper(1:end-1));

for i = 1:n_y
	loghyper_temp = loghyper([1:n_nonmeanhyper, ((n_nonmeanhyper+1+(i-1)*2*order_h):(n_nonmeanhyper+i*order_h*2))]);
% 	[nlml1, dnlml1] = GPt_likelihood_gaussmean(x, t, varargin{i}, cov_funs, mean_fun, loghyper_temp);
	[nlml1] = GPt_likelihood_gaussmean(x, t, varargin{i}, cov_funs, mean_fun, loghyper_temp, n_hyper);
	nlml = nlml + nlml1;
% 	dnlml = dnlml + dnlml1;
end