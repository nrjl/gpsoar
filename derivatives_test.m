%% Derivatives test

cov_funs = {@cressie_ns1, @d_cressie_ns1};
hyper = [0.8 15 150 1];

% Check all directions, simple 1D perturbations
delta = 1e-8;

K_1 = cov_funs{1}(x, t, x, t, hyper);

dK_numeric = zeros(1, numel(hyper));
dK_function = zeros(1, numel(hyper));

for ii = 1:numel(hyper)
	hyper2 = hyper; hyper2(ii) = hyper2(ii)+delta;
	dK_numeric(ii) = (cov_funs{1}(x,t,x,t,hyper2)-K_1)/delta;
	dK_function(ii) = cov_funs{2}(x,t,ii,hyper);	
end


%% Checkface mcmagic
[nlml, dnlml] = GPt_likelihood(X_train, t_train, W_train(:,1), cov_funs, hyper);

dnlml_numeric = zeros(numel(hyper),1);

for ii = 1:numel(hyper)
	hyper2 = hyper; hyper2(ii) = hyper2(ii)+delta;
	[nlml2, dnlml2] = GPt_likelihood(X_train, t_train, W_train(:,1), cov_funs, hyper2);
	dnlml_numeric(ii) = (nlml2-nlml)/delta;
end

disp([dnlml'; dnlml_numeric']);
