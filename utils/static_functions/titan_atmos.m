function [T, p, rho, mu] = titan_atmos(h)
global zt1999
% Note: h in km, modelled for < 40 km

% T0 = 94;		% Surface temperature
% MM = (.951*28.02 + .049*16.042)/1000;	% kg/mol
% R = 8.3145;
% RR = R/MM;
% 
% T = T0-1.15*h;
% 
% 
% A = 0.72065;
% B = -0.0128873;
% C = -0.0003245;
% D = 2.50104e-6;
% E = 6.43518e-8;
% 
% rho = 10.^(A + B*h + C*h.^2 + D*h.^3 + E*h.^4);
% 
% p = rho.*RR.*T;

[full] = interp1(zt1999(:,1), zt1999(:,2:end), h);
T = full(:,7);
p = full(:,6).*1e-4;
rho = full(:,5).*1e3;


mu = 1.718e-5+(5.1e-8*(T-273));