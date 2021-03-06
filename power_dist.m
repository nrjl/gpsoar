function C = power_dist(a, b, n)

b_full = permute(repmat(b, [1, 1, size(a, 1)]), [3, 1, 2]);
a_full = permute(repmat(a, [1, 1, size(b, 1)]), [1, 3, 2]);

C = sum((abs(b_full - a_full)).^n, 3);