function [wls] = otkrig_wls(x, t, y, vario, hyper, K_inv)
%--------------------------------------------------------------------------
%
% FUNCTION:		otkrig_wls
%
% PURPOSE:		Calculate ordinary temporal kriging weighted least squares
%               
% SYNTAX:		[wls] = otkrig_wls(x, y, vario, hyper)
%
% INPUTS:		x		- input training points [n×dim]
%				y		- output training points [n×1]
%				vario	- handle to vario estimator
%
% OUTPUTS:		nlml	- negative log marginal likelihood [1×1]
%				dnlml	- derivative of nlml wrt hyperparams [n_hyper×1]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2009
%
% MODIFIED:     November 2009
%
% See also:		otkrig_predict, ols
%--------------------------------------------------------------------------

n_x = size(x, 1);

OLS = ols_est(y);

RSS = 
WRSS1 = 