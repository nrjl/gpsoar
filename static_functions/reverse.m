function y = reverse(x)
%--------------------------------------------------------------------------
%
% FUNCTION:		reverse
%
% PURPOSE:		reverse a vector              
%              
% SYNTAX:		y = reverse(x)
%				y = reverse(x, dim)
%
% INPUTS:		x   - input matrix
%				dim - dimension to reverse along
%
% OUTPUTS:		y   - output matrix
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2007
%
% See also:		
%--------------------------------------------------------------------------

y = x;

for i = 1:length(x)
    y(i) = x(length(x)+1-i);
end    