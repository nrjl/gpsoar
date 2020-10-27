function fitness = fitness_uncertainty(pos, V, att, edot, lookahead, target_pos, E0)

global m g GR_approx V_stall

start_pos = pos(:,1,1,1);
current_pos = squeeze(pos(:,:,:,end));

Efinal = m*g*squeeze(-pos(3,:,:,end) + start_pos(3)) + 0.5*m*V(:,:,end).^2 - E0;		% [ntf,ntf] matrix
Efinal = Efinal - exp(V_stall-V(:,:,end));
a_vector = (target_pos - start_pos);
a_len = norm(a_vector);

p_size = size(current_pos);
b = squeeze(current_pos - repmat(start_pos, [1, p_size(2:3)]));	% [3,ntf,ntf]

track_distance = squeeze(dot(b, repmat(a_vector/a_len, [1, p_size(2:3)]),1));	% [ntf,ntf] matrix

E_required = m*g*(squeeze(current_pos(3,:,:) - target_pos(3)) + ...
	norm(a_vector(1:2))/GR_approx);

E_weight = (E_required >= 0)*.7 + (E_required < 0)*.4;

fitness = E_weight.*(Efinal + 0.5*Edot*lookahead) + (1-E_weight).*m*g*track_distance/GR_approx;