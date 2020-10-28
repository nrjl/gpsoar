function [K] = nick_covariance(x1, x2, hyper)

sigma_f = hyper(1);
len = hyper(2);
dist = square_dist(x1, x2);

K = sigma_f^2*(1 - 0.5*dist.^2/(len^4) + ...
	5/6*(quad_dist(x1, x2)/(len^4)) - ...
	4/3*(sqrt(dist)/len));

K = K.*(sqrt(dist) < len);