function full = segmenter(lin_sec, sections)
%--------------------------------------------------------------------------
%
% FUNCTION:		segmenter
%
% PURPOSE:		Cut trapezoidal sections into segments
%               
% INPUTS:		lin_sec	- original quarter chord points [4 x n] (x;y;z;c)
%				sections- number of segments from each of the original
%               		  sections [1 x (n-1)]
%
% OUTPUTS:		full	- full quarter chord points [4 x sum(sections)+1]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		August 2007
%
% MODIFIED:     September 2007
%
% See also:		blade_element, strip_forces
%--------------------------------------------------------------------------

% Check that there is the correct number of segment definitions
if length(sections) ~= size(lin_sec, 2)-1
    error('blade_element:segmenter:nsections', ...
        'Incorrect number of wing sections');
end


full = zeros(4, sum(sections) + 1);		% Build matrix of full points

for i = 1:length(sections)
    nsec = sections(i);					% Number of sections
    istart = 1+sum(sections(1:i-1));	% Starting index
    istop = sum(sections(1:i));			% Stop index
    full(:, istart:istop) = lin_sec(:,i)*ones(1, nsec) + ...
        (lin_sec(:,i+1)-lin_sec(:,i))*(0:1/nsec:(nsec-1)/nsec);
end
           
full(:,end) = lin_sec(:,end);			% Cap with final point