function force_plot(poa, F, F_scale)

F_end = poa + F*F_scale;

for i = 1:size(F_end, 2)
	plot3([poa(1,i), F_end(1,i)], [poa(2,i), F_end(2,i)], ...
		[poa(3,i), F_end(3,i)], 'b--');
end