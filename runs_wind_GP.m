N_RUNS = 100;
THETA_A = zeros(3, N_RUNS);
THETA_B = zeros(3, N_RUNS);
THETA_AB = zeros(3, N_RUNS);

RMS = zeros(4, N_RUNS);
LML = zeros(3, N_RUNS);
TA = zeros(1, N_RUNS);
TB = zeros(1, N_RUNS);

for NN = 1:N_RUNS
	wind_GP;
	THETA_A(:,NN) = fmina(:);
	THETA_B(:,NN) = fminb(:);
	THETA_AB(:,NN) = fminab(:);
	RMS(:,NN) = [rms_a; rms_b; rms_aob; rms_ab];
	LML(:,NN) = [nlmla; nlmlb; nlmlab];
	TA(NN) = tta;
	TB(NN) = ttb;
end
	
%%
fprintf('\n----- GP Training results -----\n');
fprintf(['Method\t\t| length\t| sigma_f\t| sigma_n\t| log(ML)\t|| ', ...
    'RMS Error\n']);
disp(repmat('-', 1, 80));
fprintf('u (alone)\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	mean(exp(THETA_A), 2),-mean(LML(1,:)), mean(RMS(1,:)));
fprintf('v (alone)\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	mean(exp(THETA_B), 2),-mean(LML(2,:)), mean(RMS(2,:)));
fprintf('uv (total)\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	[0,0,0],-mean(LML(1,:)+LML(2,:)), mean(RMS(3,:)));
fprintf('Grouped\t\t| %3.5f\t| %3.5f\t| %3.5f\t| %3.5f\t|| %3.5f\n', ...
	mean(exp(THETA_AB), 2),-mean(LML(3,:)), mean(RMS(4,:)));
