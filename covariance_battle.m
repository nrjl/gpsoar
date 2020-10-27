% Covariance battle arena

%% Simple circle, linear time

% WITH DRIFTINESS
W = [0, 0.0];

theta = linspace(0, 4*pi, 101)';
t = linspace(0, 2, numel(theta))';

x_pos = 1*cos(theta) + t*W(1);
y_pos = 1*sin(theta) + t*W(2);

X = [x_pos(:), y_pos(:)];

figure(1); clf; scatter(X(:,1), X(:,2), 5, t);
axis equal;


%% Square distance matrix
X_dist = square_dist(X, X); t_dist = square_dist(t, t);
figure(2); clf; set(2, 'position', [50, 500, 1000, 400]);
subplot(1,2,1); imagesc(max(X_dist(:)) - X_dist); axis equal; axis off; title('Spatial distance, {\it ||h||}')
subplot(1,2,2); imagesc(max(t_dist(:)) - t_dist); axis equal; axis off; title('Temporal distance, {\it |u|}')

%% Covariance
% cov_fun = @cressie_ns3;
% hyper = [1 1 1]; %[sf, lt, l];

cov_fun = @cressie_ns5;
hyper = [1, 1, 1, 5];
% hyper = [1 50*1.41 100*1.41 10];

% cov_fun =@square_expt;
% hyper = log([0.9, 0.7, 1]);

% cov_fun = @cressie_ns1;
% hyper = [1 1 1];

% Covariance variation as a function of spatial and temporal lag
h = linspace(0, 2, 100)';
u = linspace(0, 2, numel(h))';

h1 = [h(1)*ones(size(h)), zeros(size(h))];
h2 = [h, zeros(size(h))];

u1 = u;
u2 = u(1)*ones(size(u));

h_dist = square_dist(h1, h2);
u_dist = square_dist(u1, u2);

KK = cov_fun(h1, u1, h2, u2, hyper, h_dist, u_dist);
figure(5); clf; subplot(1,2,1);
[cc, hh] = contourf(h, u, KK); clabel(cc, hh);
axis tight;  daspect([1,1,1]);
xlabel('||h||'); ylabel('|u|');


K = cov_fun(X, t, X, t, hyper, X_dist, t_dist);
subplot(1,2,2); %set(gcf, 'position', [1050, 500, 500, 400]);
imagesc(K); axis equal; axis off;

%% Drifting covariance
X_distdrift = square_distdrift(X, t, X, t, W);
h_distdrift = square_distdrift(h1, u1, h2, u2, W);

KK = cov_fun(h1, u1, h2, u2, hyper, h_distdrift, u_dist);
figure(6); subplot(1,2,1);
[cc, hh] = contourf(h, u, KK); clabel(cc, hh);
axis tight;  daspect([1,1,1]);
xlabel('||h||'); ylabel('|u|');

K = cov_fun(X, t, X, t, hyper, X_distdrift, t_dist);
subplot(1,2,2);
imagesc(K); axis equal; axis off;


%% Seperable combination
hh = [log([1, 1, 1, 1]), 0, W];

KK = square_exptdrift(X, t, X, t, hh);
figure(7); clf;
imagesc(KK);




