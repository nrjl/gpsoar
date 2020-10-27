function [M] = spin_video(h, npoints)

set(h, 'Position', [100 100 640 480]);
[AZ, EL] = view;
axis vis3d

AZ_range = rectify(linspace(AZ, AZ+360, npoints+1), 'deg');

for i = 1:length(AZ_range)
	view(AZ_range(i), EL);
	drawnow;
	M(i) = getframe(gcf);
end
	

