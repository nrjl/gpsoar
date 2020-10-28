function [pos_full, att_full, V_full] = wind_sim3_control(start_pos, start_att, V0, dt, controls, lookahead, W_handle, ntf, turbulence, t0)
%--------------------------------------------------------------------------
%
% FUNCTION:		wind_sim3_control
%
% PURPOSE:		Solve 6DOF gliding equations of motion for a limited
%				control input
%
% SYNTAX:		[pos, att, V] = wind_sim3(start_pos, start_att, V0, dt, 
%					ntf, lookahead, W_handle)
%				
%				[pos, att, V, edot] = wind_sim3(...)
%
% INPUTS:		start_pos	- [3�1] Initial position (x;y;z)
%				start_att	- [3�1] Initial attitude (phi, gamma_a, psi_a)
%				V0			- [1�1] Initial airspeed
%				dt			- [1�1] Integration timestep
%				ntf			- [1�1] Control number (number of bank angles
%								and climb angles tested 
%				lookahead	- [1�1] Simulation time
%				W_handle	- Function handle to wind function of the form:
%								[W, Jw] = W_handle(pos)
%
% OUTPUTS:		pos			- [3�ntf�ntf�nt] Matrix of position during sim
%				att			- [3�ntf�ntf�nt] Matrix of attitude during sim
%				V			- [ntf�ntf�nt]   Matrix of velocity during sim
%				edot		- [ntf�ntf] Matrix of terminal energy rates
%
%				NOTE: First ntf dimension (rows) is roll, second dimension
%						(columns) is pitch
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		December 2008
%
% MODIFIED:     January 2009
%
% See also:		run_wind_sim3, wind_sim
%--------------------------------------------------------------------------
global g
global Cd0 S AR e m Nmax Nmin CL_max dphi_dt_max phi_max gamma_max
% global pitch_coeff

% ----- Timing ----- %
nt = floor(lookahead/dt);
n_segments = size(controls, 2);

% ----- Check for turbulence ----- %
% Turbulence must be in form of [u v w p q r], 6�n;
if nargin == 8
	turbulence = zeros(6, nt*n_segments);
	t0 = 0;
elseif nargin == 9
	t0 = 0;
end

pos_full = zeros(3, nt*n_segments);
att_full = zeros(3, nt*n_segments);
V_full = zeros(1, nt*n_segments);

current_pos = start_pos;
current_att = start_att;
current_V = V0;

for p = 1:n_segments

	pos = zeros(3, nt);
	pos(:,1) = current_pos;

	V = zeros(1, nt);
	V(1) = current_V;

	att = zeros(3, nt);
	att(:,1) = current_att;
	
	i = controls(1,p);
	j = controls(2,p);
	
% 	tfrac_pitch = (linspace(-1,1,ntf)-polyval(pitch_coeff, current_V))*0.5/lookahead;
	tfrac_pitch = (linspace(-1,1,ntf)-current_att(2)/(gamma_max))*0.5/lookahead;
	kstop_pitch = fix(abs(tfrac_pitch)*(nt-1));
	
	tfrac_roll  = (linspace(-1,1,ntf)-current_att(1)/(phi_max))*0.5/lookahead;
	kstop_roll  = fix(abs(tfrac_roll)*(nt-1));


	for k = 1:nt-1
		dphi_dt = (k < kstop_roll(i))*dphi_dt_max*sign(tfrac_roll(i));

% 		[T, P, rho] = titan_atmos(-pos(3,k)./1000);
		[T, P, rho] = atmos(-pos(3,k));
		V_a     = V(k);
		phi_a   = att(1,k);
		gamma_a = att(2,k);
		psi_a   = att(3,k);

% 		rho = (pos(3,i,j,k)>0)*1030 + (pos(3,i,j,k)<=0)*rho0;
		current_t = t0+(p-1)*lookahead + (k-1)*dt;
		qS = .5*rho*V_a*V_a*S;
		[W, Jw] = W_handle(pos(:,k), current_t);

		xdot = V_a*cos(gamma_a)*cos(psi_a) + W(1) + turbulence(1, ((p-1)*nt+k));
		ydot = V_a*cos(gamma_a)*sin(psi_a) + W(2) + turbulence(2, ((p-1)*nt+k));
		zdot = -V_a*sin(gamma_a) + W(3) + turbulence(3, ((p-1)*nt+k));
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
			elseif L/(m*g) > Nmax
				CL_a = Nmax*m*g/qS;
			elseif L/(m*g) < Nmin
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

		V(k+1) = V(k) + dV_dt*dt;
		pos(:,k+1) = pos(:,k) + Rdot*dt;
		att(:,k+1) = att(:,k) + dt*([dphi_dt; dgamma_dt; dpsi_dt] + ...
			calc_Ceb(att(:,k))'*turbulence(4:6, ((p-1)*nt+k)));
	end
	
	pos_full(:,(p-1)*nt+1:p*nt) = pos;
	att_full(:,(p-1)*nt+1:p*nt) = att;
	V_full((p-1)*nt+1:p*nt) = V;
	
	current_pos = pos(:,end);
	current_att = att(:,end);
	current_V = V(:,end);
end