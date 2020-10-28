load s2048_profile.mat
%%
top = s2048(1:81,:);
bottom = s2048(81:160,:);

x = linspace(0, 1, 80);
x = x(:);
x(1) = top(81,1);
ytop = interp1(top(:,1), top(:,2), x);
ybottom = interp1(bottom(:,1), bottom(:,2), x);

meanline = (ytop+ybottom)/2;

%%
figure(5); clf; hold on;
plot(top(:,1), top(:,2), 'b-', bottom(:,1), bottom(:,2), 'r-');
plot(x, meanline, 'g-');
axis equal;

%%
n_lep = 5;

phi_le = atan(mean((meanline(2:n_lep) - meanline(1))./(x(2:n_lep)-x(1))));

n_tep = 5;
phi_te = abs(atan(mean((meanline(end-n_tep:end-1) - meanline(end))./...
	(x(end-n_tep:end-1) - x(end)))));

r_le = 0.5/100;
theta = linspace(0, 2*pi, 180);
xc = r_le + cos(theta)*r_le;
yc = sin(theta)*r_le;
plot(xc, yc, 'y');