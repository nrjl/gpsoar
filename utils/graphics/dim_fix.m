function X = dim_fix(X)

X = diag([1 1 1 1 1 1 1 -1 -1 1 -1 -1])*X;
X(7,:) = rectify(X(7,:) + pi);