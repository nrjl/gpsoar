function slip_handle = slip_indicator(X, beta_max)

if any(X(1:3))
	uvw = X(1:3) - calc_Ceb(X(7:9))'*wind_field(X(10:12));
	beta = asind(uvw(2) / sqrt(uvw(1)^2 + uvw(2)^2 + uvw(3)^2));
else
	beta = 0;
end

beta_inrange = (beta_max - abs(beta))>0.5;
beta = beta_inrange*beta + (beta_max-0.5)*sign(beta)*~beta_inrange;

if beta_inrange
	col = 'k';
else
	col = 'r';
end

slip_handle = fill([0.5 0 -0.5 0 0.5]+beta, [0 0.5 0 -0.5 0], col);