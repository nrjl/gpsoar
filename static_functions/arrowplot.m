function h_out = arrowplot(x, y, z)
% Plot arrow with midpoint arrowhead

% Unscaled arrow along x axis
a_prime = [0, .5, .40, .40, .5;
			0, 0, .1, -.1, 0;
			0, 0, 0, 0, 0];

a_full = zeros(3, 5*(numel(x)-1)+1);

for ii = 1:numel(x)-1
	dx = x(ii+1) - x(ii);
	dy = y(ii+1) - y(ii);
	dz = z(ii+1) - z(ii);
	
	rr = sqrt(dx^2 + dy^2 + dz^2);
	
	if (~dx && ~dy)
		psi = 0;
	else
		psi = atan2(dy, dx);
	end
	
	theta = -asin(dz/rr);
		
	a_full(:,(5*(ii-1)+1):(5*ii)) = [x(ii); y(ii); z(ii)]*[1 1 1 1 1] + ...
		calc_Ceb(0, theta, psi)*(diag([rr 1 1])*a_prime);
	
end

a_full(:,end) = [x(end); y(end); z(end)];

h_out = plot3(a_full(1,:), a_full(2,:), a_full(3,:));