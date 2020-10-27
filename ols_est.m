function gamma_hat = ols_est(y)
%--------------------------------------------------------------------------
%
% FUNCTION:		ols_est
%
% PURPOSE:		Ordinary least squares variogram estimator
%               
% SYNTAX:		gamma_hat = ols_est(y)
%
% INPUTS:		y		- output training points [n×1]
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
% See also:		otkrig_wls
%--------------------------------------------------------------------------
s_dist = square_dist(y, y);
gamma_hat = 1/(2*numel(s_dist))*sum(s_dist(:));