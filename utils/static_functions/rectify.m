function out = rectify(in, units)
%--------------------------------------------------------------------------
%
% FUNCTION:		rectify
%
% PURPOSE:		Rectify an angle into domain (-pi, pi]
%               
% SYNTAX:		out = rectify(in, units)
%
% INPUTS:		in      - input angle
%               units   - units specification ('deg' for degrees, default
%                         is radians)
%
% OUTPUTS:		out     - output angle in specified units
%
% AUTHOR:		Nicholas Lawrance
%
% DATE:			April 2007
%
% See also:		
%--------------------------------------------------------------------------

lim = pi;                   % Default to radians

if nargin == 2              % Check whether degrees were specified
    if strcmp(units, 'deg')
        lim = 180;
    end
end

temp = in/(2*lim);          % n full rotations

out = in - fix(temp + sign(temp)/2)*2*lim;

if out/lim == -1            % Check limit values (-pi, pi]
    out = -1*out;
end