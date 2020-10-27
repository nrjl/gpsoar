function [full_controls] = calculate_path_GP(start_pos, start_att, ...
		V0, target_pos, W_handle, J_handle, tf, lookahead, dt, ntf, t0)
% Output full_controls is 2xlookahead (row 1 roll, row 2 pitch)
global g m max_climb V_stall
colours = [1 0.6 0.6; 0.6 1 0.6; 0.6 0.6 1; 1 0.6 1; 1 1 0.6; 0.6 1 1];
if nargin ~= 11
	t0 = 0;
end

fitness_f = @fitness_autoexplore;

% --- Initialisation variables --- %
E0 = m*g*-start_pos(3) + 0.5*m*V0*V0;

% IMPORTANT -  ntf is the number of roll and pitch angles to attempt
% IMPORTANT -  n_select is the number of path segments selected for next
% round
n_select = 3;	

% n_sections is the total number of sections in the full flight
n_sections = tf/lookahead;

% Handles to each section in the figure and the previous section used
h_full = zeros(n_sections, n_select);
previous_sector = zeros(n_sections, n_select);
optimum_control = zeros(2, n_sections, n_select);

% First calculation of path options. Will return ntf*ntf resulting segments
[pos1, att1, V1, edot1, R_explore1] = wind_sim3GP(start_pos, start_att, V0, dt, ntf, ...
	lookahead, W_handle, J_handle, t0);

% Rank these results based on the fitness function
E1 = m*g*squeeze(-pos1(3,:,:,end)) + 0.5*m*V1(:,:,end).^2 - E0;		% [ntf,ntf] matrix
% E1 = E1 - 3*exp(V_stall-V1(:,:,end));
E1 = E1 - 2*exp(V_stall+2 - V1(:,:,end)) - 2*exp(V1(:,:,end)-16);
% Was 22 m/s max
% Potential energy reward in field uncertainty

fitness1 = fitness_f(E1, edot1, lookahead, start_pos, pos1(:,:,:,end), target_pos, R_explore1); 	% [ntf*ntf,1] vector
[fitness_sorted, fitness_ranking] = sort(fitness1(:), 1, 'descend');

% searh matrices have dimensions:
% 1 - [data dimensions (3 for 3D data; (x,y,z), (phi,theta,psi), etc)]
% 2 - [ntf*n_select normally roll but in this case roll repeated by the
% target number of selections
pos_search = zeros(size(pos1).*[1 n_select 1 1]);
att_search = zeros(size(att1).*[1 n_select 1 1]);
V_search   = zeros(size(V1).*[n_select 1 1]);
edot_search= zeros(ntf*n_select, ntf);
fitness_search = zeros(ntf*n_select, ntf);

goal_reached = 0;

% pos_full = zeros(3, n_select, size(pos1, 4), size(pos1, 4));
% att_full = pos_full;
% V_full   = zeros(n_select, ntf, size(pos1, 4));
% edot_full= zeros(n_select, ntf, size(pos1, 4));

% Run for total number of sections, nt is the current section number (0 is
% the initial run from above)
for nt = 1:n_sections-1
	
	% Need to evaluate if the next best segment is valid, otherwise keep
	% moving through the list until n_select segments are found
	n_found = 0; i = 1;
	while n_found < n_select
		
		% Climb limit check
		while i < ntf*ntf*n_select && ~goal_reached
			[j,k] = ind2sub([ntf*(1 + (nt>1)*(n_select-1)),ntf], fitness_ranking(i));
			if abs(att1(2,j,k,end)) > max_climb
				i = i+1;
			elseif V1(j,k,end) < 3
				i = i+1;
% 			elseif 	sqrt(sum((pos1(:,j,k,end) - target_pos).^2)) < lookahead*V1(j,k,end)
% 				goal_reached = 1; break;
			else
				n_found = n_found+1;
				break
			end
		end
% 		if goal_reached; break; end;
		
		% At this point, i is the currently selected segment, n_found is
		% how many have already been found (including the current one) and
		% j,k are the indices identifying roll and pitch numbers
		% respectively.
		
% 		pos_full(:,n_found,:,nt) = pos1(:, j, k, :);
% 		att_full(:,n_found,:,nt) = att1(:, j, k, :);
% 		V_full(n_found, :, nt) = V1(j,k,:);
% 		edot_full(n_found, :, nt) = edot1(j,k,:);
		
		previous_sector(nt, n_found) = floor((j-1)/ntf)+1;
		optimum_control(:, nt, n_found) = [rem(j-1,ntf)+1;k];
		
% 		h_full(nt, n_found) = ...
% 			plot3(squeeze(pos1(1,j,k,:)), squeeze(pos1(2,j,k,:)), ...
% 			squeeze(pos1(3,j,k,:)), 'Color', colours(n_found,:));
% 		drawnow; %pause(0.1);
		
		% Extrapolate from the current segment. ntf*ntf possibilities are
		% found for each target segment from the previous round.
		t0 = t0+lookahead;
		
		[pos2, att2, V2, edot2, R_explore2] = wind_sim3GP(pos1(:,j,k,end), ...
			att1(:,j,k,end), V1(j,k,end), dt, ntf, lookahead, W_handle, J_handle, t0);
		istart = (n_found-1)*ntf+1;
		istop  = n_found*ntf;
		
		% The results are stored in the XXX_search matrices along the
		% second dimension. There are ntf*ntf*n_select total segments in
		% each round. In this stage the [ntf,ntf] results from the current
		% path extension are stacked on top of each other into a
		% [ndim, ntf*n_select, ntf] matrix (ndim is the number of data
		% dimensions)
		pos_search(:,istart:istop,:,:) = pos2;
		att_search(:,istart:istop,:,:) = att2;
		V_search(istart:istop,:,:) = V2;
		edot_search(istart:istop,:) = edot2;
		
% 		E2 = m*g*squeeze(pos1(3,j,k,end)-pos2(3,:,:,end)) + ...
% 			0.5*m*(squeeze(V2(:,:,end)).^2-V1(j,k,end).^2);	% [ntf,ntf]
		E2 = m/nt*(g*squeeze(start_pos(3)-pos2(3,:,:,end)) + ...
 			0.5*(squeeze(V2(:,:,end)).^2-V0.^2));	% [ntf,ntf]
		E2 = E2 - 2*exp(V_stall+2-V2(:,:,end)) - 2*exp(V2(:,:,end)-22);
		
		fitness_search(istart:istop,:) = fitness_f(E2, edot2, lookahead,...
			pos1(:,j,k,end), squeeze(pos2(:,:,:,end)), target_pos, R_explore2);
		i = i+1;
	end
	
	
	% Once all possibilites are found, rank the candidates using the
	% fitness function
	[fitness_sorted, fitness_ranking] = sort(fitness_search(:), 1, 'descend');
	pos1 = pos_search;
	att1 = att_search;
	V1 = V_search;
% 	if goal_reached; break; end;
end

%% Plot last paths
nt = nt + 1;
for i = 1:n_select
	[j,k] = ind2sub([ntf*(1 + (nt>1)*(n_select-1)),ntf], fitness_ranking(i));
% 	h_full(nt, i) = plot3(squeeze(pos1(1,j,k,:)), squeeze(pos1(2,j,k,:)), ...
% 		squeeze(pos1(3,j,k,:)), 'Color', colours(i,:));
	previous_sector(nt, i) = floor((j-1)/ntf)+1;
	optimum_control(:, nt, i) = [rem(j-1,ntf)+1;k];
end

%%
% final_energy = m*g*squeeze(-pos_search(3,:,:,end)) + ...
% 	0.5*m*(squeeze(V_search(:,:,end)).^2);	% [ntf,ntf]
current_sector = 1;
full_controls = zeros(2, nt);

for i = nt:-1:1
% 	set(h_full(i, current_sector), 'Color', 'b', 'LineWidth', 1);
	full_controls(:,i) = optimum_control(:, i, current_sector);
	current_sector = previous_sector(i,current_sector);
end

% set(h_full(1, current_sector), 'Color', 'r', 'LineWidth', 2);
%%
% 
% axis equal; axis tight;
% 
% xlabel('X'); ylabel('Y'); zlabel('Z');
% view(3); set(gca, 'Zdir', 'reverse');