function K = neural_net_cov(x1, x2, hyper, varargin)
% Neural net covariance function implementation
% 
% K(x,z) = sigf.^2*asin((2*x'*Sig*z)/sqrt((1+2*x'*Sig*x)*(1+2*z'*Sig*z))
% 
% Hyperparameters: [sigf, std0, std1, std2, ... , stdD]
% Where sigf is overall magnitude modifier (squared to variance), stdk
% is the standard deviation (squared to variance) of each input dimension
if nargin == 5
	x2 = hyper;
	hyper = varargin{2};
end


n_x1 = size(x1, 1);
n_x2 = size(x2, 1);

x1_tilde = [ones(n_x1, 1), x1];
x2_tilde = [ones(n_x2, 1), x2];

Sigma = diag(hyper(2:end).^2, 0);

K = (hyper(1).^2)*asin(2*x1_tilde*Sigma*x2_tilde'./...
	sqrt((1+2*sum(x1_tilde*Sigma.*x1_tilde, 2))*(1+2*sum(x2_tilde*Sigma.*x2_tilde, 2)')));