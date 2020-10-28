%% Wind field crap
% V = 2;
% lam = 100;
% K = -0.17873;
% K = K*2;
% Wx		= 0;
% Wy		= 0;
% dWx_dx	= 0;
% dWx_dy	= 0;
% dWy_dx	= 0;
% dWy_dy	= 0;
% 	V_wg_e(1,:) = Wx + dWx_dx*locs(1,:) + dWx_dy*-locs(3,:);
% 	V_wg_e(3,:) = Wy + dWy_dx*locs(1,:) + dWy_dy*-locs(3,:);


% 	alt = -(100+locs(3,:));
	
% 	V_wg_e(1,:) = (alt*-0.75.*(alt>0).*(alt < 20)...
% 		- 15*(alt >= 20));
% 
% 	V_wg_e(1,:) = alt*-0.75.*(alt>0);

% 	V_wg_e(1,:) = (alt>2).*-15;
	
% 	V_wg_e(1,:) = -pohlhausen(alt, 15, 20, 0);
% 	V_wg_e(1,:) = -log_profile(alt, 15, 20, 0.1)-ones(1,size(locs,2));
% 	V_wg_e(1:2,:) = [cos(psi_w); sin(psi_w)]*V_wg_e(1,:);

% 	V_wg_e(1,:) = -2*ones(size(locs(1,:)));
% 	V_wg_e(1,:) = -V*cos(atan(pi/5*cos(2*pi/lam*locs(1,:))));
% 	V_wg_e(2,:) = V*sin(atan(pi/5*sin(2*pi/lam*locs(1,:))));

% 	V_wg_e(1,:) = sqrt(2*locs(1,:)*K + 12^2) - 12;
% 	V_wg_e(3,:) = -0.218629222513356*ones(1, size(locs, 2));
% 	V_wg_e(3,:) = -1*ones(1, size(locs, 2));

% 	for i = 1:size(locs, 2)
% 		if (locs(1, i) > 30) && (locs(1, i) < 70)
% 			V_wg_e(2,i) = -3;
% 		elseif (locs(1, i) > 100) && (locs(1, i) < 140)
% 			V_wg_e(2,i) = 5;
% 		end
% 	end

% 	V_wg_e = torus_thermal(locs, 3, 100, 2, [200; 0; -100]);