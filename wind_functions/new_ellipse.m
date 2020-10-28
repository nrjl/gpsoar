% new ellipse

a = 2;

x = linspace(0, 1, 20);

z = a*sqrt(1 - x.^2);

dzdx = -a*x./(sqrt(1 - x.^2));

t = atan2(z,(a*x));
tprime = atan2(z,(a^2*x));

vzvx = -1./(tan(tprime));

x2 = cos(t);
z2 = a*sin(t);

figure(6); clf; hold on;
plot(x, z, 'b-', x2, z2, 'r.');
axis equal;

u1 = cos(dzdx);
v1 = sin(dzdx);

u2 = sin(tprime);
v2 = -cos(tprime);

% quiver(x, z, u1, v1);
quiver(x, z, u2, v2, 0.5);