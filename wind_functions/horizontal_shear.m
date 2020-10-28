function [W, Jw] = horizontal_shear(pos, shear, limity)

W = zeros(size(pos));
W(1,:) = (pos(2,:) < limity).*(pos(2,:) > -limity).*(pos(2,:)*shear) + ...
	(pos(2,:) >= limity).*shear*limity + (pos(2,:) <= -limity).*-shear*limity;

Jw = zeros(3,3,size(pos, 2));
Jw(1,2,:) = (pos(2,:) < limity).*(pos(2,:) > -limity).*shear;