function varargout = shift_transform(X, t, p, fh_shift, shift_param, fh_in)
%--------------------------------------------------------------------------
%
% FUNCTION:		shift_transform
%
% PURPOSE:		This function transforms an input space using a supplied 
%				transformation with respect to an additional parameter and
%				returns the function handle to the original function
%               
% SYNTAX:		C = square_dist(a, b)
%				C = square_dist(a, b, l);
%
% INPUTS:		a	- first set of points [n×d]
%				b	- second set of points [k×d]
%				l	- length scale weightings [1×d] (optional)
%
% OUTPUTS:		C	- square distances [n×k]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2009
%
% MODIFIED:     December 2009
%
% See also:		
%--------------------------------------------------------------------------

X_dash = fh_shift(X, t, shift_param);
[varargout{1}, varargout{2}] = fh_in(X_dash, p);