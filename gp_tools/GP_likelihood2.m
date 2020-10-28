function [nlml, dnlml] = GP_likelihood2(x, y1, y2, cov_funs, loghyper)
%--------------------------------------------------------------------------
%		----- OBSOLETE: REPLACE WITH GP_likelihoodn -----
%
% FUNCTION:		GP_likelihood2
%
% PURPOSE:		Calculate GP negative log marginal likelihood and 
%				derivative for two inputs
%               
% SYNTAX:		[nlml, dnlml] = GP_likelihood2(x, y1, y2, cov_funs, loghyper)
%
% INPUTS:		x		- input training points [n×dim]
%				y1		- output training points 1 [n×1]
%				y2		- output training points 2 [n×1]
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
% MODIFIED:     July 2009
%
% See also:		GP_likelihood, GP_predict
%--------------------------------------------------------------------------
disp('Warning: This function is obsolete, replace with GP_likelihoodn');

[nlml1, dnlml1] = GP_likelihood(x, y1, cov_funs, loghyper);
[nlml2, dnlml2] = GP_likelihood(x, y2, cov_funs, loghyper);

nlml = nlml1 + nlml2;
dnlml = dnlml1 + dnlml2;