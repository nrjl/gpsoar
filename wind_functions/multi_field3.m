function [W, Jw] = multi_field3(X, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		multi_field3
%
% PURPOSE:		Add multiple field component functions to yield the overall
%				field with a single function
%
% SYNTAX:		[W, Jw] = multi_field3(X, f1, p1, f2, p2, ... fn, pn)
%
% INPUTS:		X	- Matrix of test points [3×n]
%				fk	- Handle to k'th wind function of the form:
%					[V, Jw] = fk(X, pk)
%				pk	- Parameter array for corresponding function fk
%
% OUTPUTS:		W	- Output wind [3×n]
%				Jw	- Output wind Jacobian [3×3×n]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2009
%
% MODIFIED:     October 2009
%
% See also:		pohlhausen3, cos_profile3, thermal_field
%--------------------------------------------------------------------------

n_components = (nargin-1)/2;
W = zeros(size(X));
Jw = zeros(3, 3, size(X, 2));

for i = 1:n_components
	[W2, Jw2] = varargin{2*i-1}(X, varargin{2*i});
	W = W + W2;
	Jw = Jw + Jw2;
end
