function heading_handle = heading_indicator(psi)

range = 30*pi/180;
divisions = 5*pi/180;
k = 1;

low_limit = ceil((psi - range/2)/divisions)*divisions;
high_limit = floor((psi + range/2)/divisions)*divisions;
heading_handle = zeros(1, range/divisions*2);

for head_line = low_limit:divisions:high_limit
	pos = (head_line - psi)/range*2;

	
	if ~mod(head_line, 2*divisions)
		heading_handle(k) = plot([pos, pos], [-1, 0.4], 'Color', [1 1 1]);
		k = k+1;
		heading_handle(k) = text(pos, 0.75, ...
			sprintf('%3.0f', rectify(head_line)*180/pi), ...
				'Color', [1 1 1], 'FontSize', 8);
		k = k+1;
	else
		heading_handle(k) = plot([pos, pos], [-1, 0], 'Color', [1 1 1]);
		k = k+1;
	end
end

heading_handle(k) = plot([0 0], [-1, 1], 'r', 'LineWidth', 2);
heading_handle = heading_handle(1:k);