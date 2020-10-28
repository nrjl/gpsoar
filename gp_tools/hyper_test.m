lt_range = [5:1:21];
l_range = [50:40:500];

L = zeros(numel(lt_range), numel(l_range));

for ii = 1:numel(lt_range)
	for jj = 1:numel(l_range)
		L(ii,jj) = GPt_likelihood(X_train, t_train, W_train(:,1), cov_funs, log([l_range(jj), lt_range(ii), 1.101, 0.054557]));
	end
end

figure(2); clf; imagesc(lt_range, l_range, L);
xlabel('lt'); ylabel('l');

[xx, yy] = meshgrid(l_range, lt_range);
figure(3); clf; surf(xx, yy, L);
xlabel('l'); ylabel('lt'); zlabel('nlml');


%%
sf_range = [0:0.5:10];
sn_range = [0.05:0.05:0.5];

L2 = zeros(numel(sf_range), numel(sn_range));

for ii = 1:numel(sf_range)
	for jj = 1:numel(sn_range)
		L2(ii,jj) = GPt_likelihood(X_train, t_train, W_train(:,1), cov_funs, [2, 50, 300, sf_range(ii), sn_range(jj)]);
	end
end

figure(4); clf; imagesc(sf_range, sn_range, L2);
xlabel('sf'); ylabel('sn');

[xx, yy] = meshgrid(sn_range, sf_range);
figure(5); clf; surf(xx, yy, L2);
xlabel('sn'); ylabel('sf'); zlabel('nlml');