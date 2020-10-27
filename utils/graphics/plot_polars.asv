%--------------------------------------------------------------------------
%
% SCRIPT:		plot_polars
%
% PURPOSE:		Plot and perform calculations for a specific wing profile
%               
% INPUTS:		Data file must contain data stored as a 3-dimensional 
%               array with 1st dimension alpha, 2nd dimension variable (see
%               list), 3rd dimension reynolds number. Reynolds numbers must
%               be stored in an a vector named Re such that length(Re) =
%               size(polars, 3). File must also contain name of profile
%               stored in a string called foilname
%
% OUTPUTS:		(to command window and graphically)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		March 2007
%
% MODIFIED:     June 2007
%
% See also:		sbxc_full, atmos
%--------------------------------------------------------------------------

%% Polars Script

%Columns:   1       2   3   4   5   6   7   8   9
%           alpha	Cl	Cd	Cm 	TU	TL	SU	SL	L/D

load S2048_control.mat %NACA001054_polars.mat

if ~(exist('polars', 'var') && exist('Re', 'var'))
    error('WingModel:plot_polars:MissingVar', ...
        'Polars matrix or Reynolds number list is missing')
end

if exist('foilname', 'var') == 1
    fprintf('\n----- %s Aerofoil Polars Analysis -----\n', foilname)
end

n_Re = size(polars, 3);
if n_Re ~= length(Re)
    error('WingModel:plot_polars:wrongdim', ...
        'Specified Reynolds numbers do not match 3rd dimension of polars')
end

n_alpha = size(polars, 1);

%% Variables

mg = 53.3958;
L = 0.290;
[T, P, rho, a, mu, nu] = atmos(100);
A = 0.9738;


%% Cl(Cd) polar

figure(1); clf; hold on;
zerodef = find(def == 0);
zerodef = zerodef(1);
%subplot(2,1,2); cla; hold on;

Cl = squeeze(polars(:,2,:,zerodef));
Cd = squeeze(polars(:,3,:,zerodef));
Cm = squeeze(polars(:,4,:,zerodef));

lines = {'b-s', 'r--o', 'g-.+', 'k:x'};

lgnd = {};
Clo = [];
Cdo = [];

min_Re = 2e5;
max_Re = 5e5;
min_Cl = min(min(Cl));
max_Cl = max(max(Cl));

for i = 1:n_Re
    plot(Cd(:,i), Cl(:,i), lines{mod(i-1, 4)+1});
    
    if (Re(i) >= min_Re) && (Re(i) <= max_Re)
        index = find(Cd(:,i) < 0.03);
        Clo = [Clo; Cl(index, i)];
        Cdo = [Cdo; Cd(index, i)];
		min_Cl = max([min_Cl, Cl(index(1), i)]);
		max_Cl = min([max_Cl, Cl(max(index), i)]);
    end
        
    lgnd{i} = sprintf('Re = %d', Re(i));
end
    
p_Cd= polyfit(Clo, Cdo, 2);
Cli = linspace(min_Cl, max_Cl, 20);
Cdi = polyval(p_Cd, Cli);
plot(Cdi, Cli, 'k-', 'linewidth', 1.5);
lgnd{i+1} = 'Best quadratic fit';

Cl_points = [min_Cl, -p_Cd(2)/2/p_Cd(1), max_Cl];
Cd_points = polyval(p_Cd, Cl_points);
plot(Cd_points, Cl_points, 'k^')

legend(lgnd, 'Location', 'best')
grid on;

ax = axis;
plot([ax(1), ax(2)], [0, 0], 'k-')
plot([0, 0], [ax(3), ax(4)], 'k-')

quadstr = sprintf(['Quadratic approximation: C_d = %0.4gC_l^2', ...
    '+ %0.4gC_l + %0.4g'], p_Cd(1), p_Cd(2), p_Cd(3));
disp(quadstr)
fprintf('Cl:\t%0.5g\t\t%0.5g\t\t%0.5g\n', Cl_points(1), Cl_points(2), Cl_points(3));
fprintf('Cd:\t%0.5g\t%0.5g\t%0.5g\n', Cd_points(1), Cd_points(2), Cd_points(3));
text(ax(1)+.03*(ax(2)-ax(1)), ax(4)-.03*(ax(4)-ax(3)), quadstr)

title(['C_{\it{l}}(C_{\it{d}}) Polar - Quadratic Fit (Re \leq', ...
    sprintf(' %g)', max_Re)])
xlabel('Drag Coefficient, C_{\it{d}}')
ylabel('Lift Coefficient, C_{\it{l}}')

%% Cl(alpha)
figure(2); clf; hold on;
%subplot(2,1,1); cla; hold on;

alfa = polars(:,1,1);

index = find((alfa >= 0).*(alfa <= 5));
a1 = alfa(index)*pi/180;
c1 = Cl(index,1);

p_a = polyfit(a1, c1, 1);

a2 = [min(alfa), max(alfa)]*pi/180;
c2 = polyval(p_a, a2);

for i = 1:n_Re
    plot(alfa, Cl(:,i), lines{mod(i-1, 4)+1});
end

plot(a2*180/pi, c2, 'k-', 'linewidth', 1.5)
lgnd{i+1} =  'Best linear fit';

legend(lgnd, 'Location', 'best')
grid on;

ax = axis;
plot([ax(1), ax(2)], [0, 0], 'k-')
plot([0, 0], [ax(3), ax(4)], 'k-')

linstr = sprintf('Linear approximation: C_l = %0.5g + %0.5g', p_a(2), p_a(1));
fprintf(1, '\n%s*alpha\n', linstr)
text(ax(1)+.03*(ax(2)-ax(1)), ax(4)-.03*(ax(4)-ax(3)), [linstr, '\alpha'])

title('Lift coefficient vs Angle of Incidence \alpha')
xlabel('Angle of Incidence, \alpha (deg)')
ylabel('Lift Coefficient, C_{\it{l}}')

%% Pitching moment
figure(3); clf; hold on;

index = find((alfa >= -2).*(alfa <= 6));
a1 = alfa(index)*pi/180;
m1 = Cm(index,1);

p_m = polyfit(a1, m1, 1);

a2 = [min(alfa), max(alfa)]*pi/180;
m2 = polyval(p_m, a2);

for i = 1:n_Re
    plot(alfa, Cm(:,i), lines{mod(i-1, 4)+1});
end

plot(a2*180/pi, m2, 'k-', 'linewidth', 1.5);

legend(lgnd, 'Location', 'best')
grid on;
ax = axis;

linstr = sprintf('Linear approximation: C_m = %0.5g + %0.5g', p_m(2), p_m(1));
fprintf(1, '\n%s*alpha\n', linstr)
text(ax(1)+.03*(ax(2)-ax(1)), ax(4)-.03*(ax(4)-ax(3)), [linstr, '\alpha'])
plot([ax(1), ax(2)], [0, 0], 'k-')

title('Pitching Moment vs Angle of Incidence \alpha')
xlabel('Angle of Incidence, \alpha (deg)')
ylabel('Pitching Moment, C_{\it{m}}')


%% L/D(alpha)
L_D = Cl./Cd;

figure(4); clf; hold on;

for i = 1:n_Re
    plot(alfa, L_D(:,i), lines{mod(i-1, 4)+1});
end

%lgnd = lgnd(1:i);

legend(lgnd{1:i}, 'Location', 'best')
grid on;

% ax = axis;
% plot([ax(1), ax(2)], [0, 0], 'k-')
% plot([0, 0], [ax(3), ax(4)], 'k-')

title('Lift-to-drag ratio vs Angle of Incidence \alpha')
xlabel('Angle of Incidence, \alpha (deg)')
ylabel('Lift to Drag Ratio, ^L/_D')


%% Cruise velocity

v2 = 11:0.1:14;

Re2 = v2*L/nu;

Cl2 = mg./(.5*rho*v2.^2*A);

alpha2 = zeros(1, length(v2));
L_D2 = alpha2;
D = alpha2;

for i = 1:length(v2)
    [y, i_Re] = min(abs(Re - Re2(i)));
    i_Re = [i_Re, i_Re+sign(Re2(i) - Re(i_Re))];
    Re_r = Re(i_Re);
    
    % Linear interplotation between Reynolds numbers
    Cl_fit = interp1(Re_r, Cl(:, i_Re)', Re2(i), 'linear');
    
    
    alpha2(i) = interp1(Cl_fit(index), alfa(index), Cl2(i), 'linear');
    
    L_D2(i) = interp2(Re, alfa, L_D, Re2(i), alpha2(i), 'linear');
end

[max_LD, i_max] = max(L_D2);
v_max = v2(i_max);
a_max = alpha2(i_max);

fprintf(1, '\nMaximum performance obtained for a cruise speed of %0.4g m/s\n', ...
    v_max);
fprintf(1, 'Requires a main wing angle of incidence of %0.4g degrees\n', a_max)
fprintf(1, 'Gives an overall L/D of %0.5g\n', max_LD);

figure(5); clf;
subplot(2,1,1); hold on;
plot(v2, alpha2, 'b-');
plot(v_max, a_max, 'ro');
text(v_max+.1, a_max, [sprintf('V = %0.4g,' , v_max), '\alpha = ',...
    sprintf('%0.4g', a_max)]);
xlabel('Cruise velocity (m/s)');
ylabel('Wing angle of incidence (deg)');
grid on;


subplot(2,1,2); hold on;
plot(v2, L_D2, 'b-');
plot(v_max, max_LD, 'ro');
text(v_max+.1, max_LD, sprintf('V = %0.4g, L/D = %0.4g', v_max, max_LD));
xlabel('Cruise velocity (m/s)');
ylabel('Lift to drag ratio L/D');
grid on;

figure(6); clf;
subplot(2,1,1);
surf(Re, alfa, Cl)
xlabel('Re')
ylabel('\alpha (deg)')
zlabel('C_L')
view(-77,22);

subplot(2,1,2);
surf(Re, alfa, L_D)
xlabel('Re')
ylabel('\alpha (deg)')
zlabel('L/D')
view(-77,22);


figure(7); clf;
subplot(2,1,1);
[C, H] = contourf(Re, alfa, L_D, 10);
set(H, 'LineStyle', 'none');
xlabel('Re')
ylabel('\alpha (deg)')
hold on;
plot(Re2, alpha2, 'k-', 'LineWidth', 2);

subplot(2,1,2);
surf(Re, alfa, L_D)
xlabel('Re')
ylabel('\alpha (deg)')
zlabel('L/D')
alpha(0.5);
hold on;
plot3(Re2, alpha2, L_D2, 'k-o', 'LineWidth', 2);


% ------------------------- %
% ------- SCRAPBOOK ------- %
% ------------------------- %
%% Cruise Velocity

% figure(4); clf;
% subplot(2,1,1); hold on;
% subplot(2,1,2); hold on;
% 
% v = Cl;
% Re_LD = Cl;
% 
% for i = 1:n_Re
%     v(:,i) = sqrt(2*mg./(rho*A*Cl(:,i)));
%     
%     % Re_LD(:,i) = v*L/nu;
%     D = 0.5*rho*A*v(:,i).^2.*Cd(:,i);
%     
%     subplot(2,1,1); hold on;
%     plot(alfa, v(:,i), lines{mod(i-1, 4)+1});
%     
% end
% 
% subplot(2,1,2); hold on;
% plot(alfa, D, lines{mod(i-1, 4)+1});
%     
% % lgnd = lgnd(1:i);
% 
% legend(lgnd{1:i}, 'Location', 'best')
% grid on;
% 
% subplot(2,1,1);
% title('Cruise speed vs Angle of Incidence \alpha')
% xlabel('Angle of Incidence, \alpha (deg)')
% ylabel('Cruise speed, m/s')
% grid on;
% 
% subplot(2,1,2);
% title('Induced Drag vs Angle of Incidence \alpha')
% xlabel('Angle of Incidence, \alpha (deg)')
% ylabel('Total Wing Induced Drag, N')
% grid on;


% Cl_fit = (Re2(i) - *(Cl([i_Cl-1, i_Cl, i_Cl+1], i_Re) + ...
%    Cl([i_Cl-1, i_Cl, i_Cl+1], i_Re +sign(Re2(i)-Re(i_Re))));
%     linp = polyfit(Cl([i_Cl-1, i_Cl, i_Cl+1], i_Re), ...
%         L_D([i_Cl-1, i_Cl, i_Cl+1], i_Re), 1);
%     
%     L_D2(i) = polyval(linp, Cl2(i));
%     
%     D(i) = mg/(L_D2(i));

%     L_D2(i) = interp1(Cl_fit, L_D, Cl2(i), 'linear');
    
%     linp = polyfit(Cl([i_Cl-1, i_Cl, i_Cl+1], i_Re), ...
%         L_D([i_Cl-1, i_Cl, i_Cl+1], i_Re), 1);
%     
%     L_D2(i) = polyval(linp, Cl2(i));
%     
%     D(i) = mg/(L_D2(i));