V_core = 1;
R = 1;
k = 3;

% Plot a toroidal thermal field
X = linspace(-2*R, 2*R, 25);
Y = linspace(-2*R, 2*R, 25);
Z = linspace(-k*R, k*R, 7);

[XX, YY, ZZ] = meshgrid(X, Y, Z);

[U,V,W] = torus_thermal(XX, YY, ZZ, V_core, R, k, [0,0,0]);

%%
figure(5); clf;
quiver3(XX, YY, ZZ, U, V, W);
axis equal;
set(gca, 'ZDir', 'rev')

%%
figure(6); clf; hold on;
wscale = (Z(2) - Z(1))/V_core/1.5;
hs = zeros(1, numel(Z));
for i = 1:numel(Z)
	hs(i) = surf(XX(:,:,1), YY(:,:,1), W(:,:,i)*wscale + Z(i), W(:,:,i));
end
set(gca, 'ZDir', 'rev')
set(hs, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
xlabel('X');
ylabel('Y');
zlabel('Z');
colorbar; view(3);

%% Awesome ring plot thing
r = 100;		% Plot radius
h = k*r;

n_t = 16;		% Number of points around each ring
t = linspace(-pi, pi, n_t);

n_u = 15;		% Number of rings
u = linspace(-pi, pi, n_u);

figure(7);
hold on;

X2 = zeros(1, n_u*n_t);
Y2 = X2; Z2 = X2;

for i = 1:length(u)
    X = (R + r/2*cos(t)).*cos(u(i));
    Y = (R + r/2*cos(t)).*sin(u(i));
    Z = h/2*sin(t);
    
    dex = (i-1)*n_t+1:i*n_t;
    X2(dex) = X;
    Y2(dex) = Y;
    Z2(dex) = Z;
% 	plot3(X, Y, Z, 'Color', [0.1, 0.5, 0.1]);
end

W2 = torus_thermal([X2;Y2;Z2], V_core, R, k, [0;0;0]);
figure(7); hold on;
quiver3(X2, Y2, Z2, W2(1,:), W2(2,:), W2(3,:)) ;
axis equal;
set(gca, 'ZDir', 'rev')