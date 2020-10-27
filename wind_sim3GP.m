function [pos, att, V, Edot, R_explore] = ...
	wind_sim3GP(start_pos, start_att, start_V, dt, ntf, lookahead, W_handle, J_handle, t0)
%--------------------------------------------------------------------------
%
% FUNCTION:		wind_sim3GP
%
% PURPOSE:		Solve 6DOF gliding equations of motion for a limited
%				control input with GP 
%
% SYNTAX:		[pos, att, V] = wind_sim3GP(start_pos, start_att, start_V, 
%					dt, ntf, lookahead, W_handle, J_handle)
%				
%				[pos, att, V, edot] = wind_sim3(...)
%
% INPUTS:		start_pos	- [3×1] Initial position (x;y;z)
%				start_att	- [3×1] Initial attitude (phi, gamma_a, psi_a)
%				start_V		- [1×1] Initial airspeed
%				dt			- [1×1] Integration timestep
%				ntf			- [1×1] Control number (number of bank angles
%								and climb angles tested)
%				lookahead	- [1×1] Simulation time
%				W_handle	- Wind function of the form: W = W_handle(pos)
%				J_handle	- Wind Jacobian:	J = J_handle(pos)
%
% OUTPUTS:		pos			- [3×ntf×ntf×nt] Matrix of position during sim
%				att			- [3×ntf×ntf×nt] Matrix of attitude during sim
%				V			- [ntf×ntf×nt]   Matrix of velocity during sim
%				edot		- [ntf×ntf] Matrix of terminal energy rates
%
%				NOTE: First ntf dimension (rows) is roll, second dimension
%						(columns) is pitch
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		December 2008
%
% MODIFIED:     September 2009
%
% See also:		run_wind_sim3, wind_sim
%--------------------------------------------------------------------------
global g
global Cd0 S AR e m Nmax Nmin CL_max dphi_dt_max phi_max gamma_max
global x_limits y_limits z_limits
global pitch_coeff

if nargin ~= 9
	t0 = 0;
end

% ----- Timing ----- %
t = 0:dt:lookahead;
nt = length(t);

tfrac_pitch = (linspace(-1,1,ntf)-start_att(2)/(gamma_max))*0.5/lookahead;
% tfrac_pitch = (linspace(-1,1,ntf)-polyval(pitch_coeff, start_V))*0.5/lookahead;
tfrac_roll  = (linspace(-1,1,ntf)-start_att(1)/(phi_max))*0.5/lookahead;
% ntf = length(tfrac);

% ----- Initialisation ----- %
% [T0, P0, rho0] = titan_atmos(-start_pos(3)./1000);
[T0, P0, rho0] = atmos(-start_pos(3));

pos = zeros(3, ntf, ntf, nt);
pos(:,:,:,1) = reshape(start_pos(:)*ones(1, ntf*ntf), [3,ntf,ntf,1]);

V = zeros(ntf, ntf, nt);
V(:,:,1) = start_V*ones(ntf,ntf,1);

att = zeros(3, ntf, ntf, nt);
att(:,:,:,1) = reshape(start_att(:)*ones(1, ntf*ntf), [3,ntf,ntf,1]);

Edot = zeros(ntf, ntf);

kstop_pitch = fix(abs(tfrac_pitch)*(nt-1));
kstop_roll  = fix(abs(tfrac_roll)*(nt-1));

for i = 1:ntf

	for j = 1:ntf		
		
		for k = 1:nt-1
			
			dphi_dt = (k < kstop_roll(i))*dphi_dt_max*sign(tfrac_roll(i));
			
			V_a     = V(i,j,k);
			phi_a   = att(1,i,j,k);
			gamma_a = att(2,i,j,k);
			psi_a   = att(3,i,j,k);

			rho = (pos(3,i,j,k)>0)*1030 + (pos(3,i,j,k)<=0)*rho0;
			qS = .5*rho*V_a*V_a*S;
			W = W_handle(pos(:,i,j,k)', t0+t(k));
			Jw = J_handle(pos(:,i,j,k)', t0+t(k));

			xdot = V_a*cos(gamma_a)*cos(psi_a) + W(1);
			ydot = V_a*cos(gamma_a)*sin(psi_a) + W(2);
			zdot = -V_a*sin(gamma_a) + W(3);
			Rdot = [xdot; ydot; zdot];

			if k >= kstop_pitch(j) % If pitch rate is given
				dgamma_dt = 0;

				L = m/cos(phi_a)*(V_a*dgamma_dt + g*cos(gamma_a) - ...
					[[cos(psi_a), sin(psi_a)]*sin(gamma_a), cos(gamma_a)]*Jw*Rdot);

				CL_a = L/qS;

				if CL_a > CL_max
					CL_a = CL_max;
				elseif CL_a < -CL_max
					CL_a = -CL_max;
				elseif L/(m*g) > Nmax;
					CL_a = Nmax*m*g/qS;
				elseif L/(m*g) < Nmin;
					CL_a = Nmin*m*g/qS;
				end
			else
				CL_a = (j < kstop_pitch(j))*((tfrac_pitch(j) < 0)*Nmin + (tfrac_pitch(j)>=0)*Nmax)*m*g/qS;
			end

			L = CL_a*qS;

			dgamma_dt = 1/V_a*(L/m*cos(phi_a) - g*cos(gamma_a) + ...
				[[cos(psi_a), sin(psi_a)]*sin(gamma_a), cos(gamma_a)]*Jw*Rdot);

			D = Cd0*qS + L*L/(pi*AR*e*qS);

			dV_dt = -D/m - g*sin(gamma_a) - ...
				[[cos(psi_a), sin(psi_a)]*cos(gamma_a), -sin(gamma_a)]*Jw*Rdot;

			dpsi_dt = 1/(V_a*cos(gamma_a))*(L/m*sin(phi_a) + ...
				[sin(psi_a), -cos(psi_a), 0]*Jw*Rdot);

			V(i,j,k+1) = V(i,j,k) + dV_dt*dt;
			pos(:,i,j,k+1) = pos(:,i,j,k) + Rdot*dt;
			att(:,i,j,k+1) = att(:,i,j,k) + [dphi_dt; dgamma_dt; dpsi_dt]*dt;
		end
		Edot(i,j) = m*(g*-zdot + V(i,j,k)*dV_dt);
	end
end

% R_explore output
if nargout == 5
	
	end_pos = reshape(pos(:,:,:,end), 3, ntf*ntf);
	end_V = reshape(V(:,:,end), 1, ntf*ntf);
	start_V = start_V*ones(1, ntf*ntf);
	
	[end_W, end_VW] = W_handle(end_pos', (t0+lookahead)*ones(ntf*ntf, 1));
	end_stdW = repmat(sqrt(end_VW'), 3, 1);
	
	[start_W, start_VW] = W_handle(start_pos', t0);
	start_stdW = repmat(sqrt(start_VW), 3, 1);
		
	R_explore = explore_reward(repmat(start_pos, 1, ntf*ntf), end_pos, ...
		start_V, end_V, repmat(start_stdW, 1, ntf*ntf), end_stdW, lookahead);
	
	low_lim = repmat([x_limits(1); y_limits(1); z_limits(1)], 1, ntf*ntf);
	high_lim= repmat([x_limits(2); y_limits(2); z_limits(2)], 1, ntf*ntf);
	R_explore = R_explore.*all((end_pos > low_lim) & (end_pos < high_lim), 1);
	R_explore = reshape(R_explore, ntf, ntf);
end

%	

% 	end_att = reshape(att(:,:,:,end), 3, ntf*ntf);
% 	
% 	for i = 1:ntf*ntf
% 		R_explore(i) = explore_reward(start_pos, start_V, start_att, ...
% 			end_pos(:,i), end_V(:,i), end_att(:,i), start_stdW, end_stdW(:,i), lookahead);
% 	end
	 

% %% Extra outputs
% 	X = reshape(pos(1,:,:,:), [ntf*ntf, nt])';
% 	Y = reshape(pos(2,:,:,:), [ntf*ntf, nt])';
% 	Z = reshape(pos(3,:,:,:), [ntf*ntf, nt])';
% 	varargout{1} = plot3(X, Y, Z, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');

% E = m*g*Y(:,end) + .5*m*(V(:,end).^2);
% E_adj = 0.5*m*g*V(:,end).*exp(2*(V_stall - V(:,end)));
% [Emax, imax] = max(E(:,end) - E_adj);
% 
% if nargout > 4
% 	varargout{1} = t;
% 	varargout{2} = tfrac;
% 	varargout{3} = L_m;
% 	varargout{4} = D_m;
% 	varargout{5} = imax;
% end
