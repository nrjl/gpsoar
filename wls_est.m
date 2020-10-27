function gamma_bar = wls_est(y)
%--------------------------------------------------------------------------
%
% FUNCTION:		wls_est
%
% PURPOSE:		Weighted least squares variogram estimator
%               
% SYNTAX:		gamma_bar = wls_est(y)
%
% INPUTS:		y		- output training points [n×1]
%
% OUTPUTS:		gamma_bar - Cressie weighted least squares
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2009
%
% MODIFIED:     November 2009
%
% See also:		ols_est
%--------------------------------------------------------------------------
s_dist = square_dist(y, y);
gamma_bar = (1/numel(s_dist)*sum(s_dist(:).^(1/4)))^4/...
			(0.914 + 0.988/numel(s_dist));