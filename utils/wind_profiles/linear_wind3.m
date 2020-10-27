function [W, Jw] = linear_wind3(pos, W0, Jw)
% [W] = linear_wind3(pos, W0, Jw)

W = W0*ones(1, size(pos,2)) + Jw*pos;
Jw = repmat(Jw, [1, 1, size(pos,2)]);