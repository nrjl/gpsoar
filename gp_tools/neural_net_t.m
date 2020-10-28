function Kxt = neural_net_t(x1, t1, x2, t2, hyper)
% Hyper: [sigmaf, sigma0, sigmax, sigmay, sigmaz, sigmat]
Kxx = neural_net_cov(x1, x2, hyper(1:5));		% Training spatial covariances
Ktt = square_exp(t1, t2, [hyper(6), 0]);		% Training temporal covariances

Kxt = Kxx.*Ktt;