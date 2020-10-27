function Kxt = seperable_cov(x1, t1, x2, t2, fun1, fun2, hyper1, hyper2)
% Kxt = seperable_cov(x1, t1, x2, t2, fun1, fun2, hyper1, hyper2)
Kxx = fun1(x1, x2, hyper1);		% Training spatial covariances
Ktt = fun2(t1, t2, hyper2);		% Training temporal covariances

Kxt = Kxx.*Ktt;

