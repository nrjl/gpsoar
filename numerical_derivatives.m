function J = numerical_derivatives(x, y, x_test, cov_funs, loghyper, Ki)


% V0 = GP_predict(x, y, x_test, cov_funs{1}, loghyper, Ki);

% Get numerical derivatives
n_test = size(x_test, 1);
delta = (max(x) - min(x))*1e-8;
zerox = zeros(size(x_test));
deltax = delta(1)*ones(n_test, 1);
deltay = delta(2)*ones(n_test, 1);
deltaz = delta(3)*ones(n_test, 1);


% dx
dx = zerox; dx(:,1) = deltax;
Vxh = GP_predict(x, y, x_test+dx, cov_funs{1}, loghyper, Ki);
Vxl = GP_predict(x, y, x_test-dx, cov_funs{1}, loghyper, Ki);
dV_dx = (Vxh - Vxl)/(2*delta(1));

% dy
dy = zerox; dy(:,2) = deltay;
Vyh = GP_predict(x, y, x_test+dy, cov_funs{1}, loghyper, Ki);
Vyl = GP_predict(x, y, x_test-dy, cov_funs{1}, loghyper, Ki);
dV_dy = (Vyh - Vyl)/(2*delta(2));

% dz
dz = zerox; dz(:,3) = deltaz;
Vzh = GP_predict(x, y, x_test+dz, cov_funs{1}, loghyper, Ki);
Vzl = GP_predict(x, y, x_test-dz, cov_funs{1}, loghyper, Ki);
dV_dz = (Vzh - Vzl)/(2*delta(3));

J_numeric = [dV_dx; dV_dy; dV_dz];

% GP Prediction
J = dGP_predict(x, y, x_test, [1,2,3], cov_funs{2}, loghyper, Ki);

% Abs Difference
disp(abs(J-J_numeric));
