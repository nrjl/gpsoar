function R_explore = explore_reward(p0, p1, v0, v1, std0_est, std1_est, delta_t)

global m g
p_vec = abs(p1 - p0);
p_len = sqrt(sum((p1-p0).^2, 1));

W_vec = [(std1_est(1:2,:) + std0_est(1:2,:)); abs(std1_est(3,:) - std0_est(3,:))];

% Scalar projection of the 
max_deltaV = dot(W_vec, p_vec, 1)./p_len;

R_explore = m*g*(std1_est(3,:) + std0_est(3,:))*delta_t./2 + ...
	.5*m*((v0+v1).*max_deltaV + max_deltaV.^2);

%% SCRAP
% a_mean = atan2(sin(a0)+sin(a1), cos(a0)+cos(a1)); % Angular mean
% gam = atan((p0(3,:)-p1(3,:))/sqrt((p1(1,:) - p0(1,:)).^2 + (p1(2,:) - p0(2,:)).^2));
% psi = atan2( p1(2,:)-p0(2,:) , p1(1,:)-p0(1,:) );
% 
% s2_gam = (sin(gam)).^2;
% c2_gam = (cos(gam)).^2;
% s2_psi = (sin(psi)).^2;
% c2_psi = (cos(psi)).^2;
% 
% dWx = std1_est(1,:) + std0_est(1,:);
% dWy = std1_est(2,:) + std0_est(2,:);
% dWz = std1_est(3,:) - std0_est(3,:);
% 
% max_deltaV = sqrt(c2_psi.*c2_gam.*dWx.^2 + s2_psi.*c2_gam.*dWy.^2 + ...
% 	s2_gam.*dWz.^2);
%
%
% R_explore = explore_reward(p0, V0, a0, p1, V1, a1, std0_est, std1_est,
% delta_t)
%
% s_gamma = sign(a_mean(2));
% 
% gradient_directions = [.5, -sign(sin(2*a_mean(3))), sign(cos(a_mean(3))).*s_gamma;
% 					 0, .5, sign(sin(a_mean(3))).*s_gamma;
% 					 0, 0, .5];
% gradient_directions = gradient_directions + gradient_directions';
% 
% dW = [(std1_est(1:2) + std0_est(1:2)); (std1_est(3) - std0_est(3))];
% Jw_magnitude = abs([dW./(p1(1)-p0(1)), dW./(p1(2)-p0(2)), dW./(p1(3)-p0(3))]);
% Jw_magnitude(isnan(Jw_magnitude) | isinf(Jw_magnitude)) = 0;
% Jw_best = gradient_directions.*Jw_magnitude;
% 
% delta_W = (Jw_best*(p1 - p0));
% 
% 
% Va = (V0+V1)./2*[cos(a_mean(3))*cos(a_mean(2)); sin(a_mean(3))*cos(a_mean(2)); -sin(a_mean(2))];
% 
% E_dot_est = m*(g*(std0_est(3) + std1_est(3))./2 + ...
% 	-Va'*Jw_best*(Va + delta_W));
% 
% R_explore = E_dot_est*delta_t;