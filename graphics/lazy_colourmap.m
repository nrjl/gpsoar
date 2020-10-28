% my lazy colourmap checker

colour = @(dex) [0, 255, 0]./255 + [1,-1,0]*dex + [0,0,-2*abs(dex-0.5)+1];

cmat = zeros(100, 3);
di = 0.01;
figure(2); clf; hold on;

for i = 0:di:1
	cmat(round(i/di)+1,:) = colour(i);
	fill([i-di, i, i, i-di], [0,0,di,di], colour(i))
end

plot(0:di:1, [cmat(:,3), cmat(:,2), cmat(:,1)]*di);