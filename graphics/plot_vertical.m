function h = plot_vertical(x, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		plot_vertical
%
% PURPOSE:		Plot vertical lines between current axis limits 
%					(mainly used for mode changes in plot_states)
%
% SYNTAX:		plot_vertical(x)	- Plots vertical lines at specified x			
%				plot_vertical(x, ls)
%				plot_vertical(x, ls, col, width)
%
% INPUTS:		x	- x positions of lines (will use all points specified)
%				ls  - line style (standard styles, '-', '-.', etc.)
%				col	- colour (standard colors 'b', 'g', etc or [1x3] RGB
%
% OUTPUTS:		(graphical - current axis between y-limits)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		May 2008
%
% MODIFIED:     July 2008
%
% See also:		plot_states
%--------------------------------------------------------------------------
ax = axis;
ls = '-';
col = [255 170 0]./255;
wid = 1;

if nargin == 2
	ls = varargin{1};
elseif nargin == 3
	ls = varargin{1};
	col = varargin{2};
elseif nargin == 4
	ls = varargin{1};
	col = varargin{2};
	wid = varargin{3};
end

x = x(:);
nline = length(x);
h = zeros(nline, 1);

for i = 1:nline
	h(i) = plot([x(i), x(i)], [ax(3), ax(4)], ls, 'Color', col);
end