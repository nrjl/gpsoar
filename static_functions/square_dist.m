function C = square_dist(a, b, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		square_dist
%
% PURPOSE:		Calculate square distances between sets of points
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

if nargin >2
	l = varargin{1};
	a = a./repmat(l, [size(a,1), 1]);
	b = b./repmat(l, [size(b,1), 1]);
end

b_full = permute(repmat(b, [1, 1, size(a, 1)]), [3, 1, 2]);
a_full = permute(repmat(a, [1, 1, size(b, 1)]), [1, 3, 2]);

C = sum((b_full - a_full).^2, 3);