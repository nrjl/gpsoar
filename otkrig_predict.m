function [f_star, V] = otkrig_predict(x, t, y, x_star, t_star, vario, hyper, K_inv)
%--------------------------------------------------------------------------
%
% FUNCTION:		otkrig_predict
%
% PURPOSE:		Make an ordinary temporal Kriging prediction and return 
%				mean and covariance at specified target points
%               
% SYNTAX:		[f_star] = otkrig_predict(x, t, y, xs, ts, vario, hyper)
%				[f_star, V] = otkrig_predict(x, t, y, xs, ts, vario, hyper)
%
%				To just provide inv(K) use empty test point set:
%				[Ki] = GP_predict(x, t, y, [], [], vario, hyper)
%
%				To solve for a known Ki = inv(K):
%				[f_star, V] = ...
%					otkrig_predict(x, t, y, xs, ts, vario, hyper, Ki)
%
% INPUTS:		x		- training point locations
%				t		- training point times
%				y		- training point output values
%				xs		- test point locations
%				ts		- test point times
%				vario	- handle to covariance function
%				hyper	- hyperparameters
%				Ki		- Inverted Ki matrix
%
% OUTPUTS:		f_star	- test point output
%				V		- point variances
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2009
%
% MODIFIED:     November 2009
%
% See also:		
%--------------------------------------------------------------------------
n_x = size(x, 1);
	
if nargin ~= 8
	
	K = vario(x, t, x, t, hyper);		% Training point covariances
	
	K_mat = ones(n_x+1, n_x+1); K_mat(n_x+1, n_x+1) = 0;
	K_mat(1:n_x, 1:n_x) = K;
	
	% Cholesky decomposition:
	% K*A = y
	% K = L*L'
	% Use Matlab Cholesky decomp to get L (lower trangular)
	% Lz = y	=>	z = L\y
	% L'x = z	=>	A = L'\z
	%			=>	A = L'\(L\y)
	
	L = chol(K_mat, 'lower');
	K_inv = (L'\(L\eye(n_x)));
end

% OUTPUT
if isempty(x_star)
	f_star = K_inv;
else
	n_x_star = size(x_star, 1);
	v_star = [vario(x, t, x_star, t_star, hyper); ones(1, n_x_star)];
	weights = K_inv*v_star;
	
	f_star = weights(1:n_x-1, :)'*y;

	if nargout == 2
		% Variances of test points
		V = weights'*v_star;
	end
end