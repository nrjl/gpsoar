function [W, Jw] = sinusoidal_wind3(pos, V, lambda)

[Wx, Wz, dWx_dx, dWz_dx] = sinusoidal_wind(pos(1,:), V, lambda);

W = [Wx; zeros(size(Wx)); Wz];

Jw = zeros(3,3,size(pos,2));
Jw(1,1,:) = dWx_dx;
Jw(3,1,:) = dWz_dx;