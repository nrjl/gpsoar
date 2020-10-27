% figure(1); clf;
% figure(2); clf;
% figure(3); clf;
% figure(4); clf;

% ellipse test

r = 100;		% Inner radius (actually minor diameter of small ellipses)
k = 3;			% Height ratio (1 = circle, >1 for ellipses, ratio of axes)
h = k*r;
R = 100;			% Big radius
vc = 3;			% Core velocity

n_t = 21;

t = linspace(-pi, pi, n_t);
x = R + r/2*cos(t);
y = h/2*sin(t);

%%
figure(11);
hold on;
% plot(x, y, '-b', x, y, 'r.');

tprime = atan2(y*r^2,(x-R)*h^2);

a = vel_set(x, vc, R, tprime);
u2 = a.*sin(tprime);
v2 = -a.*cos(tprime);

quiver(x, y, u2, v2, 0.5);
axis equal;


figure(12);
plot(x, u2, '-b', x, v2, 'r:', x, a, 'g--');
legend('x-vel', 'z-vel', 'total vel');

%%
n_u = 20;
u = linspace(-pi, pi, n_u);

figure(13);
hold on;

X2 = zeros(1, n_u*n_t);
Y2 = X2; Z2 = X2;
U2 = X2; V2 = X2; W2 = X2;


for i = 1:length(u)
    X = (R + r/2*cos(t)).*cos(u(i));
    Y = (R + r/2*cos(t)).*sin(u(i));
    Z = h/2*sin(t);
    
    dex = (i-1)*n_t+1:i*n_t;
    X2(dex) = X;
    Y2(dex) = Y;
    Z2(dex) = Z;
    
    tprime = atan2(Z*r^2, (sqrt(X.^2 + Y.^2)-R)*h^2);
    
    a = vel_set(sqrt(X.^2 + Y.^2), vc, R, tprime);
    U2(dex) = a.*sin(tprime).*cos(u(i));
    V2(dex) = a.*sin(tprime).*sin(u(i));
    W2(dex) = -a.*cos(tprime);

    plot3(X, Y, Z, 'g-');
end

axis equal

figure(14);
hold on;
quiver3(X2, Y2, Z2, U2, V2, W2);
axis equal;


%% Distribution field

n_p = 7;

tr = linspace(-2*R, 2*R, n_p*n_p);
dex = 1:n_p:n_p*n_p;

[X3, Y3, Z3] = meshgrid(tr, tr, tr*k/2);
d = sqrt(X3.^2 + Y3.^2);

rall = 2/k*sqrt((d - R).^2*k^2  + Z3.^2);

u2 = atan2(Y3,X3);
t2 = atan2(Z3,(d-R)*k^2).*(X3 | Y3) + pi*~(X3 | Y3);

% a = (2 - sqrt(X3.^2 + Y3.^2 + Z3.^2));
% a = a.*(a > 0);
% a = ones(size(X3));

a = vel_set(d, vc, R, t2);

U3 = a.*sin(t2).*cos(u2);
V3 = a.*sin(t2).*sin(u2);
W3 = -a.*cos(t2);

figure(15); clf;
quiver3(X3(dex, dex, dex), Y3(dex, dex, dex), Z3(dex, dex, dex), ...
    U3(dex, dex, dex), V3(dex, dex, dex), W3(dex, dex, dex), 1);
axis equal;

figure(16); clf; hold on;

Wfull = zeros(n_p*n_p);
Wflat = W3(:,:,1);
Index = find(sqrt(X3(:,:,1).^2 + Y3(:,:,1).^2) <= 2*R);

for i = 1:length(Index)
    Wfull(Index(i)) = Wflat(Index(i));
end

surf(X3(:,:,1), Y3(:,:,1), Wfull./vc);
view(3);
xlabel('x')
ylabel('y')
zlabel('w/w_{core}')
