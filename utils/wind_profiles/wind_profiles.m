% Wind profiling code

z = 0:0.1:20;
u_inf = 20;
delta = 10;

%% Standard profiles from boundary layer theory

zd = (z<=delta).*z;

% Logarithmic
h0 = 0.01;
u_log = log(z/h0)/log(delta/h0);

% ERF
u_erf = erf(z/delta);

% Seventh
u_sev = (zd/delta).^(1/7) + (z>delta);

% Polhausen's quartic
lambda = 0;
eta = zd/delta;
u_pol = 2*eta - 2*eta.^3 + eta.^4 + lambda/6*eta.*(1-eta).^3 + (z>delta);

% Standard quadratic profile
u_quad = (-(zd/delta-1).^2 + 1) + (z>delta);

% Sine
u_sin = sin(zd/delta*pi/2) + (z>delta);

%% Plotting
figure(3); clf;
plot(u_log, z/delta, u_erf, z/delta, u_sev, z/delta, u_pol, z/delta, u_quad, ...
    z/delta, u_sin, z/delta);
xlabel('Tangential velocity u/u_\infty');
ylabel('Height (y/\delta)');
legend('Logarithmic', 'Error function', '1/7^{th} profile', ['Polhausens',...
    ' quartic (\Lambda = ', sprintf('%1.0f)', lambda)] ,'Quadratic', 'Sine', ...
    'location', 'best');
