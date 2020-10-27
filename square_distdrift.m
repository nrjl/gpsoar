function C = square_distdrift(x1, t1, x2, t2, W)
%--------------------------------------------------------------------------
%
% FUNCTION:		square_distdrift
%
% PURPOSE:		Calculate square distances between drifted sets of points
%               
% SYNTAX:		C = square_distdrift(x1, t1, x2, t2, W)
%				C = square_distdrift(x1, t1, x2, t2, W, l)
%
% INPUTS:		x1	- first set of points [n×d]
%				t1	- first set of times [n×1]
%				x2	- second set of points [k×d]
%				t2	- second set of times [k×1]
%				l	- length scale weightings [1×d] (optional)
%
% OUTPUTS:		C	- square distances [n×k]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		September 2010
%
% MODIFIED:     September 2010
%
% See also:		
%--------------------------------------------------------------------------

n = size(x1,1); k = size(x2,1);
d = size(x1,2);

x1_full = permute(repmat(x1, [1, 1, k]), [1, 3, 2]);
x2_full = permute(repmat(x2, [1, 1, n]), [3, 1, 2]);

t_mat = repmat(t1', [k, 1]) - repmat(t2, [1, n]);
W_full = repmat(t_mat', [1,1,d]);

for ii = 1:d
	W_full(:,:,ii) = W_full(:,:,ii).*W(ii);
end

C = sum((x1_full + W_full - x2_full).^2, 3);