function plot_states(X, U, tx, tu, m, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		plot_states
%
% PURPOSE:		Plot states
%
% SYNTAX:		plot_states(X, U, tx, tu, m)
%				plot_states(X, U, tx, tu, m, oplot, col, ls, modetimes)
%
% INPUTS:		X	- [n×m] State vector - n states, m points in time
%				U	- [j×k] Control vector - j controls, k points in time
%				tx	- [1×m] Timing vector for X
%				tu	- [1×k] Timing vector for U
%				m	- aircraft mass
%				oplot - overplot flag (1 for overplot, 0 for new plots)
%				col - plot colour (must be a string)
%				ls	- line style (string)
%
% OUTPUTS:		(graphical - figures 2-8)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2007
%
% MODIFIED:     November 2007
%
% See also:		test, control_test
%--------------------------------------------------------------------------

global g ctime
fsize = 8;

modetimes = 0; t0 = tx(1); tf = tx(end);
main_line_colour = [255 100 0  ]./255;
mode_line_colour = [255 230 100]./255;

switch nargin
	case 5
		oplot = 0;
		col = 'b';
		ls = '-';
	case 6
		oplot = varargin{1};
		cols = 'bgrmy'; col = cols(oplot+1);
		lss = {'-', '--', '-.', ':', '-'}; ls = lss{oplot+1};
	case 7
		oplot = varargin{1};
		col = varargin{2};
		ls = '-';
	case 8
		oplot = varargin{1};
		col = varargin{2};
		ls = varargin{3};
	case 9
		oplot = varargin{1};
		col = varargin{2};
		ls = varargin{3};
		modetimes = varargin{4};
end

set(0,'DefaultFigureWindowStyle','docked')

% Wind relative speeds
uvw_wind = zeros(3, size(X,2));		% Wind-relative velocity in body frame
v_earth = zeros(3, size(X,2));		% Inertial velocity in earth frame
wind_full = wind_field(X(10:12,:));

for i = 1:size(X,2)
	ctime = tx(i);
	v_earth(:,i) = calc_Ceb(X(7:9,i))*X(1:3,i);
	uvw_wind(:,i) = X(1:3,i) - calc_Ceb(X(7:9,i))'*wind_full(:,i);
end

% Body speeds
figure(2);
set(gcf, 'Name', 'Body Air Speeds')
subplot(3,1,1);
if (~oplot); cla; end; hold on;
plot(tx, uvw_wind(1,:), [col, ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('forward body speed, u (m/s)', 'FontSize', fsize); grid on;

subplot(3,1,2);
if (~oplot); cla; end; hold on;
plot(tx, uvw_wind(2,:), [col, ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('side body speed, v (m/s)', 'FontSize', fsize); grid on;

subplot(3,1,3);
if (~oplot); cla; end; hold on;
plot(tx, uvw_wind(3,:), [col, ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('vertical body speed, w (m/s)', 'FontSize', fsize); grid on;

% Body rates and angles
figure(3);
set(gcf, 'Name', 'Body rates and angles')
subplot(2,1,1);
if (~oplot); cla; end; hold on;
plot(tx, X(4,:)*180/pi, ['b', ls], ...
	 tx, X(5,:)*180/pi, ['r', ls], ...
	 tx, X(6,:)*180/pi, ['g', ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('rotation rates (deg/s)', 'FontSize', fsize); grid on;
h_l = legend('Roll Rate, p', 'Pitch Rate, q', 'Yaw Rate, r');
set(h_l, 'location', 'best', 'FontSize', fsize);

subplot(2,1,2);
if (~oplot); cla; end; hold on;
plot(tx, X(7,:)*180/pi, ['b', ls], ...
	 tx, X(8,:)*180/pi, ['r', ls], ...
	 tx, X(9,:)*180/pi, ['g', ls], ...
	 tx, atan2(v_earth(2,:), v_earth(1,:))*180/pi, ['m', ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('Euler angles (deg)', 'FontSize', fsize); grid on;
h_l = legend('Bank Angle \phi', 'Pitch Angle, \theta', ...
	'Heading Angle, \psi', 'Bearing');
set(h_l, 'location', 'best', 'FontSize', fsize);

% Aero conditions
figure(4); set(gcf, 'Name', 'Airspeed and Aero angles');
V = sqrt(sum(uvw_wind(1:3,:).^2, 1));
alfa = atan2(uvw_wind(3,:), uvw_wind(1,:))*180/pi;
beta = asind(uvw_wind(2,:)./V);

subplot(2,1,1);
if (~oplot); cla; end; hold on;
plot(tx, V, [col, ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('Airspeed (m/s)', 'FontSize', fsize); grid on;

subplot(2,1,2);
if (~oplot); cla; end; hold on;
plot(tx, alfa, ['b', ls], tx, beta, ['r', ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('Aero angle (deg)', 'FontSize', fsize); grid on;
h_l = legend('Angle of attack, \alpha', 'Sideslip, \beta');
set(h_l, 'location', 'best', 'FontSize', fsize);

% Energy
figure(5); set(gcf, 'Name', 'Energy');

subplot(2,1,1);
if (~oplot); cla; end; hold on;
Ek = 0.5*m*V.^2;
Eks = Ek - Ek(1);
Ep = m*g*-X(12,:);
Eps = Ep - Ep(1);
Et = Eks + Eps;
plot(tx, Eks, ['b', ls], tx, Eps, ['r', ls], tx, Et, ['g', ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('Energy Change, \DeltaE (J)', 'FontSize', fsize); grid on;
h_l = legend('Kinetic Energy', 'Potential Energy', 'Total Energy');
set(h_l, 'location', 'best', 'FontSize', fsize);

subplot(2,1,2);
if (~oplot); cla; end; hold on;

diff_t = diff(tx);
plot(tx(1:end-1) + diff_t/2, diff(Eks)./diff_t, ['b', ls], ...
	 tx(1:end-1) + diff_t/2, diff(Eps)./diff_t, ['r', ls], ...
	 tx(1:end-1) + diff_t/2, diff(Et)./diff_t, ['g', ls]);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize); 
ylabel('Energy rate, Edot (J/s, W)', 'FontSize', fsize); grid on;
axis tight; a_lim = axis;
a_lim = a_lim.*[1,1,1.2,1.2] + [0, 0, -1, 1]*oplot*0.1*(a_lim(4)-a_lim(3));
axis(a_lim);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
text(a_lim(1)+0.1*(a_lim(2)-a_lim(1)), a_lim(3)+0.1*(a_lim(4)-a_lim(3)),...
	sprintf('Average energy rate = %0.5g J/s', (Et(end)-Et(1))/tx(end)),...
	'FontSize', fsize);

% Position Plot
figure(6);
set(gcf, 'Name', 'Position');
subplot(3,1,1);
if (~oplot); cla; end; hold on;

plot(tx, X(10,:), [col, ls]);

set(gca, 'Xlim', [t0, tf]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
xlabel('Time (s)', 'FontSize', fsize); 
ylabel('X Position (m)', 'FontSize', fsize); grid on;

subplot(3,1,2);
if (~oplot); cla; end; hold on;

plot(tx, X(11,:), [col, ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('Y Position (m)', 'FontSize', fsize); grid on;

subplot(3,1,3);
if (~oplot); cla; end; hold on;

plot(tx, -X(12,:), [col, ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('Altitude (m)', 'FontSize', fsize); grid on;

% Control Deflections
figure(7);

U2 = reshape([U(:,1:end-1); U(:,1:end-1)], size(U, 1), 2*size(U, 2)-2);
tu2 = reshape([tu(1:end-1); tu(2:end) - diff(tu)/1000], 1, 2*size(U, 2)-2);

if (~oplot); clf; end; hold on;
set(gcf, 'Name', 'Control deflections');

plot(tu2, U2, ls);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('Control Deflection [-1,1]', 'FontSize', fsize); grid on;
h_l = legend('Elevator', 'Aileron', 'Rudder', 'Flaps');
set(h_l, 'Location', 'best', 'FontSize', fsize);

% Energy trade
figure(8);
set(gcf, 'Name', 'Energy trade');
if (~oplot); clf; end; hold on;
plot(Eks, Eps, [col, ls]); %Eks, Eps, [col,'x']
for k = 1:length(tx)
	if abs(tx(k) - 0.1*round(tx(k)/0.1)) < 1e-8
		plot(Eks(k), Eps(k), [col, '.'])
	end
end
xlabel('Kinetic Energy (J)', 'FontSize', fsize);
ylabel('Potential Energy (J)', 'FontSize', fsize); grid on;
axis equal;
ax = axis; hold on;
plot([0,  min([ax(2), -ax(3)])], [0, -min([ax(2), -ax(3)])], 'k-', ...
	 [0, -min([-ax(1), ax(4)])], [0,  min([-ax(1), ax(4)])], 'k-')

 for i = 1:numel(modetimes)
	 ind = find(tx == modetimes(i));
	 plot(Eks(ind), Eps(ind), 'ro');
 end

figure(10);
set(gcf, 'Name', 'Climb angle');
if (~oplot); clf; end; hold on;

plot(tx, asin(-v_earth(3,:)./sqrt(sum(v_earth(1:3,:).^2, 1)))*180/pi, ['r', ls]);
plot(tx, X(8,:)*180/pi - alfa, ['b', ls]);
plot_vertical(modetimes(:,2:end), ls, mode_line_colour);
plot_vertical(modetimes(:,1), ls, main_line_colour);
set(gca, 'Xlim', [t0, tf]);
xlabel('Time (s)', 'FontSize', fsize);
ylabel('Climb Angle (deg)', 'FontSize', fsize); grid on;
legend('Inertial climb angle, \gamma_i', 'Air-relative climb angle, \gamma_a')

set(0,'DefaultFigureWindowStyle','normal')