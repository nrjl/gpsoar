% 3D Dipole potential field

q = 2;
a = 1;

x = -2:0.5:2;

[X, Y, Z] = meshgrid(x, x, x);


U = q/(4*pi)*(-X.*(X.^2 + Y.^2 + (Z+a).^2).^(-3/2) + X.*(X.^2 + Y.^2 + (Z-a).^2).^(-3/2));

V = q/(4*pi)*(-Y.*(X.^2 + Y.^2 + (Z+a).^2).^(-3/2) + Y.*(X.^2 + Y.^2 + (Z-a).^2).^(-3/2));

W = q/(4*pi)*(-(Z+a).*(X.^2 + Y.^2 + (Z+a).^2).^(-3/2) + (Z-a).*(X.^2 + Y.^2 + (Z-a).^2).^(-3/2));

%%
figure(1); clf;
quiver3(X, Y, Z, U, V, W);
grid off;

%%
eps = 1;

U2 = 3*X.*Z.*eps./(4*pi*(X.^2 + Y.^2 + Z.^2).^(5/2));

V2 = 3*Y.*Z.*eps./(4*pi*(X.^2 + Y.^2 + Z.^2).^(5/2));

W2 = -eps.*((X.^2 + Y.^2 + Z.^2)-3*Z.^2)./(4*pi*(X.^2 + Y.^2 + Z.^2).^(5/2));

%%
figure(2); clf;
quiver3(X, Y, Z, U2, V2, W2, 1.5);
grid off;