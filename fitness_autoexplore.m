function fitness = fitness_autoexplore(Efinal, Edot, lookahead, start_pos, current_pos, target_pos, R_explore)

global m g GR_approx
global K_explore

% a is vector from start pos to target pos
% b is vector from current pos to target pos

a_vector = (target_pos-start_pos);
a_len = norm(a_vector);

p_size = size(current_pos);

% b = squeeze(current_pos - repmat(start_pos, [1, p_size(2:3)]));	% [3,ntf,ntf] vectors from start to current
% track_distance = squeeze(dot(b, repmat(a_vector/a_len, [1, p_size(2:3)]),1));	% [ntf,ntf] matrix

b = squeeze(repmat(target_pos, [1, p_size(2:3)]) - current_pos); 	% [3,ntf,ntf] vectors from target to current
b_len = squeeze(sqrt(sum(b.^2, 1)));
track_distance = (a_len - b_len);

% E_required = m*g*(squeeze(current_pos(3,:,:) - target_pos(3)) + ...
% 	norm(a_vector(1:2))/GR_approx);

min_E_weight = 0.2;
E_weight_100 = 0.5;

weight_auto = (1-min_E_weight)*exp(log((E_weight_100-min_E_weight)/(1-min_E_weight))/100*b_len) + min_E_weight;

% E_weight = (E_required >= 0).*weight_auto + (E_required < 0)*.4;
E_weight = weight_auto;

fitness = E_weight.*(0.5*Efinal + 0.5*Edot*lookahead) + (1-E_weight).*m*g*track_distance/GR_approx + K_explore*R_explore;