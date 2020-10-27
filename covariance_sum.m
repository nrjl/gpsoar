function [K] = covariance_sum(x1, t1, x2, t2, hyper, fun_handles, n_hyper)
% [K] = covariance_sum(x1, t1, x2, t2, hyper, fun_handles, n_hyper)
% n_hyper is a list of the number of hyperparameters of the corresponding
% covariance function. Any arguments for the mean function must be appended
% AFTER the other covariance functions, but should still have a 
% corresponding entry in the n_hyper vector, but will not have a 
% fun_handles entry (i.e numel(fun_handles) = numel(h_hyper) or 
% numel(fun_handles) = numel(h_hyper)-1.
%
% Example:  For a sum of two covariance functions with 4 HPs each, and a 
% mean function with 6 HPs the corresponding n_hyper is [4 4 6].
n_functions = numel(fun_handles);

K = zeros(size(x1, 1), size(x2, 1));

h_start = 1;

for i = 1:n_functions
	fh_covariance = fun_handles{i};
	h_end = (h_start+n_hyper(i)-1);
	Kt = fh_covariance(x1, t1, x2, t2, hyper(h_start:h_end));
	K = K + Kt;
	
	h_start = h_end+1;
end