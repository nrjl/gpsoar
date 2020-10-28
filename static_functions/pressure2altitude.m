function z = pressure2altitude(P)
%--------------------------------------------------------------------------
%
% FUNCTION:		pressure2altitude
%
% PURPOSE:		Matlab code to calculate the equivalent altitude in the 
%				International Standard Atmosphere for a given pressure. 
%				Available for altitudes 0 - 86000m.
%
%	NOTE: Need corrections to deal with places below sea level (Amsterdam)
%
%               
% SYNTAX:		z = pressure2altitude(P)
%
% INPUTS:		P - pressure (Pa) (must be a vector; column or row)
%
% OUTPUTS:		z - geometric altitude (m)
%
% AUTHOR:		Nicholas Lawrance
%
% DATE:			July 2011
%
% See also:		atmos
%--------------------------------------------------------------------------

% Constants
To = 288.15;    % K
R = 8.31432;    % N.m/mol.K
M = 0.02896442; % kg/mol (air)

E = 6356766;    % Radius of Earth (m)
g0 = 9.80665;   % Gravitational acceleration at SL
GM_R = g0*M/R;  % GM/R term

% Get Lapse rate dependence on altitude (note geopotential altitude)
h_layers = [0, 11, 20, 32, 47, 51, 71, 84.852]*1e3; % Layer heights in m
L_layers = [-6.5, 0, 1.0, 2.8, 0, -2.8, -2.0]/1e3;  % Layer lapse rates K/m

% Layer boundary temperatures (K)
T_layers = To + [0, cumsum(L_layers.*diff(h_layers))];	

% Note that this is a bit of a hack to reduce requirement for full
% calculation of layer pressures (i.e limited to 7 SF precision)
P_layers = [101325, 22632.04, 5474.878, 868.0158, 110.9058, 66.93853, ...
	3.956393, 0.37338105];

h = zeros(size(P));

for ii = 1:numel(P)
	 
	layer = sum(P(ii) <= P_layers);
	
	if layer < 1 || layer > 7
		warning('Pressure of %0.5g Pa outside ISA range', P(ii));
		h(ii) = 0;
		continue;
	end
	
	L = L_layers(layer);
	
	if L == 0
		h(ii) = h_layers(layer) - T_layers(layer)/GM_R*log(P(ii)/P_layers(layer));
		
	else
		T = T_layers(layer)*(P(ii)/P_layers(layer))^(-L/GM_R);
		h(ii) = h_layers(layer) + (T-T_layers(layer))/L;
	end
end

% Convert geopotential altitude to geometric altitude.
z = h.*E./(E-h);
	

