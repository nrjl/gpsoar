function [W, Jw] = pohlhausen3(pos, z0, u_inf, delta, lambda)

W = zeros(size(pos));
[W(1,:), dWx_dz] = pohlhausen(-pos(3,:), z0, u_inf, delta, lambda);

Jw = zeros(3,3,size(pos, 2));
% dWx/dz
Jw(1,3,:) = -dWx_dz;