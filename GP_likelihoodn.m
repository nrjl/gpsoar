function [nlml, dnlml] = GP_likelihoodn(x, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		GP_likelihoodn
%
% PURPOSE:		Calculate GP negative log marginal likelihood and 
%				derivative for n inputs
%               
% SYNTAX:		[nlml, dnlml] = ...
%					GP_likelihoodn(x, y1, y2, ... yn, cov_funs, loghyper)
%
% INPUTS:		x		- input training points [i×dim]
%				yk		- output training points k [n×1]
%				cov_funs- handles to covariance function and its derivative
%							in structure form:	{cov_fun, dcov_fun} 
%				loghyper- log of hyperparameters [n_hyper×1]
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

if nargin < 4
	error('GPNick:GP_likelihoodn:n_inputs', 'Too few input arguments');
else
	cov_funs = varargin{end-1};
	loghyper = varargin{end};
	n_y = nargin - 3;
end

nlml = 0;
dnlml = 0;

for i = 1:n_y
	[nlml1, dnlml1] = GP_likelihood(x, varargin{i}, cov_funs, loghyper);
	nlml = nlml + nlml1;
	dnlml = dnlml + dnlml1;
end