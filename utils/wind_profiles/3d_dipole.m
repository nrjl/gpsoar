% 3D Dipole potential field

q = 1;
a = 1;

x = -5:0.1:5;

[X, Y, Z] = meshgrid(x, x, x);


U = q/(4*pi)*(-X.*(X.^2 + Y.^2 + (Z+a).^2)^(-3/2) + X.*(X.^2 + Y.^2 + (Z-a).^2)^(-3/2));

V = q/(4*pi)*(-Y.*(X.^2 + Y.^2 + (Z+a).^2)^(-3/2) + X.*(X.^2 + Y.^2 + (Z-a).^2)^(-3/2));

W = q/(4*pi)*(-(Z+a).*(X.^2 + Y.^2 + (Z+a).^2)^(-3/2) + (Z-a).*(X.^2 + Y.^2 + (Z-a).^2)^(-3/2));

figure(1); clf;
quiver3(X, Y, Z, U, V, W);
