function [UVW, PQR] = dryden_TF(h, V, b, u_20, n, dt)
%--------------------------------------------------------------------------
%
% FUNCTION:		dryden_TF
%
% PURPOSE:		Generate continuous turbulence estimate from Dryden PSD for
%				low altitude (<1000 ft) flight.
%               
% SYNTAX:		[UVW, PQR] = dryden_TF(h, V, b, u_20, n, dt)
%
% INPUTS:		h		- altitude (m)
%				V		- airspeed (m/s)
%				b		- wing span (m)
%				u_20	- mean wind speed at 20ft (m/s)
%							Light turbulence: 7.617 m/s
%							Moderate turbulence: 15.234 m/s
%							Severe turbulence: 30.468 m/s
%				n		- number of time steps
%				dt		- time step increment (s)
%
% OUTPUTS:		UVW		- matrix of turbulence magnitudes (m/s) [n×3]
%				PQR		- matrix of rotational turbulence (rad/s) [n×3]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2010
%
% MODIFIED:     July 2010
%
% See also:		
%--------------------------------------------------------------------------

% Units to Imperial (ft, lb, s)
u_20 = u_20*3.281; 
h = h*3.281;
b = b*3.281;
V = V*3.281;

% Low altitude model

% Turbulence scale lengths
L_u = h/(0.177+0.000823*h).^(1.2);
L_v = L_u;
L_w = h;

% Turbulence intensities
sigma_w = 0.1*u_20;
sigma_u = sigma_w/(0.177+0.000823*h).^(0.4);
sigma_v = sigma_u;

% Dryden turbulence TFs
s = tf('s');

dryden_uTF = sigma_u*sqrt(2*L_u/(pi*V))/(1+s*L_u/V);

dryden_pTF = sigma_w*sqrt(0.8/V)*(pi/(4*b))^(1/6)/...
	(L_w^(1/3)*(1+s*4*b/pi/V));

dryden_vTF = sigma_v*sqrt(L_v/pi/V)*(1 + s*sqrt(3)*L_v/V)/(1+s*L_v/V)^2;

dryden_rTF = s/V/(1+s*3*b/(pi*V))*dryden_vTF;

dryden_wTF = sigma_w*sqrt(L_w/pi/V)*(1 + s*sqrt(3)*L_w/V)/(1+s*L_w/V)^2;

dryden_qTF = s/V/(1+s*4*b/(pi*V))*dryden_wTF;


% Convert to state space
[num, den] = tfdata(dryden_uTF, 'v');
[Au,Bu,Cu,Du] = tf2ss(num,den);

[num, den] = tfdata(dryden_vTF, 'v');
[Av,Bv,Cv,Dv] = tf2ss(num,den);

[num, den] = tfdata(dryden_wTF, 'v');
[Aw,Bw,Cw,Dw] = tf2ss(num,den);

[num, den] = tfdata(dryden_pTF, 'v');
[Ap,Bp,Cp,Dp] = tf2ss(num,den);

[num, den] = tfdata(dryden_qTF, 'v');
[Aq,Bq,Cq,Dq] = tf2ss(num,den);

[num, den] = tfdata(dryden_rTF, 'v');
[Ar,Br,Cr,Dr] = tf2ss(num,den);

% Generate time series
rr = pi*randn([6, n]);

% These are the internal state representations, need to be stored but not
% human useful; they are transferred through C and D to get actual gust
% values
xxu = zeros(size(Au,1), 1);
xxv = zeros(size(Av,1), 1);
xxw = zeros(size(Aw,1), 1);
xxp = zeros(size(Ap,1), 1);
xxq = zeros(size(Aq,1), 1);
xxr = zeros(size(Ar,1), 1);

% Preallocate gusts
UVW = zeros(n, 3);
PQR = zeros(n, 3);

% Numerical integration
for ii = 1:n
	% Newton lazy integration
	xxu = xxu + (Au*xxu + Bu*rr(1,ii))*dt;
	UVW(ii,1) = Cu*xxu;
	
	xxv = xxv + (Av*xxv + Bv*rr(2,ii))*dt;
	UVW(ii,2) = Cv*xxv;
	
	xxw = xxw + (Aw*xxw + Bw*rr(3,ii))*dt;
	UVW(ii,3) = Cw*xxw;
	
	xxp = xxp + (Ap*xxp + Bp*rr(4,ii))*dt;
	PQR(ii,1) = Cp*xxp;
	
	xxq = xxq + (Aq*xxq + Bq*rr(5,ii))*dt;
	PQR(ii,2) = Cq*xxq;
	
	xxr = xxr + (Ar*xxr + Br*rr(6,ii))*dt;
	PQR(ii,3) = Cr*xxr;
end
	
UVW = UVW./3.281;

% tt = 0:dt:((n-1)*dt);
% figure(1); plot(tt', UVW)
% figure(2); plot(tt', PQR)