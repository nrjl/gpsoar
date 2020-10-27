% Pull data out of save files
% clear all

% Get list of filenames
filelist = what('Data\JGCD_Data\WeakThermal\');
folder = filelist.path;
filenames = filelist.mat;

% Remove non save files (files without 'Save' in them)
cutlist=zeros(1, 20); cutn = 1;
for ij = 1:numel(filenames)
	if isempty(strfind(filenames{ij}, 'Save'))
		cutlist(cutn) = ij;
		cutn = cutn+1;
	end
end
filenames(cutlist(1:cutn-1)) = '';

nfiles = numel(filenames);
load([folder, filenames{1}]);

n_frames = size(t_train_full,2);
all_variance = zeros(nfiles, n_frames+1);
all_maxlift = zeros(1, nfiles);
all_winderror = zeros(nfiles, n_frames+1);

for ij = 1:nfiles
	
	load([folder, filenames{ij}]);
	
	all_variance(ij,:) = uncertainty;
	all_maxlift(ij) = max(thermals(:,1));
	
	WW = W_actual([x_grid(:), y_grid(:), z_grid(:)]', zeros(1, numel(x_grid)));
	
	for jj = 1:n_frames
		% Plot start features
		CTrain = squeeze(X_train_full(:,:,jj));
		WTrain = squeeze(W_train_full(:,:,jj));
		cc = max_obs; while CTrain(cc,1)==0; cc=cc-1; end;
		CTrain =CTrain(1:cc, :); WTrain =WTrain(1:cc, :);
		
		[WW_grid, VW_grid] = GP_predict(CTrain, WTrain, ...
			[x_grid(:), y_grid(:), z_grid(:)], cov_funs{1}, loghyper);
		all_winderror(ij,jj) = sqrt(mean(sum((WW_grid-WW').^2, 2)));
	end
	
	
	
	fprintf(1,'%s complete\n', filenames{ij});
	
end
all_winderror(:,n_frames+1) = all_winderror(:,n_frames);
all_variance(:,n_frames+1) = all_variance(:,n_frames);


%%
figure(1); clf;
dskip = 11;
t_variance = [t0:replan:(n_replan-1)*replan, tf];

boxplot(all_variance(:,1:dskip:end)./numel(VW_grid));
xlabel('Time (s)');
ylabel('Mean field variance,  \it{mean}\rm(\sigma_i^2), (m^2/s^2)');

%
figure(2); clf;
h_var = plot(t_variance, all_variance./numel(VW_grid), 'k-', 'color', .5*[1 1 1]);
h_vax = gca; hold on;
xlabel('Time (s)');
ylabel('Mean field variance,  \it{mean}\rm(\sigma_i^2), (m^2/s^2)');
% grid on;

% Plot winderror
figure(3); clf;
boxplot(all_winderror(:,1:dskip:end));
xlabel('Time (s)');
ylabel('RMS wind error  (m/s)');


figure(4); clf;
h_err = plot(t_variance, all_winderror, '-', 'color', .5*[1 1 1]);
h_eax = gca; hold on;
xlabel('Time (s)');
ylabel('RMS wind error  (m/s)');
% grid on;

% Find worst cases
[vmax, vjj] = max(all_variance(:,end));
[emax, ejj] = max(all_winderror(:,end));

% set(h_var(vjj), 'color', .2*[1 1 1], 'linewidth', 1.0);
set(h_var(ejj), 'color', .2*[1 1 1], 'linewidth', 1.0, 'linestyle', '--');
tt1 = h_var(1); h_var(1) = h_var(ejj); h_var(ejj) = tt1;
set(h_vax, 'children', h_var);

% set(h_err(vjj), 'color', .2*[1 1 1], 'linewidth', 1.0);
set(h_err(ejj), 'color', .2*[1 1 1], 'linewidth', 1.0, 'linestyle', '--');
tt1 = h_err(1); h_err(1) = h_err(ejj); h_err(ejj) = tt1;
set(h_eax, 'children', h_err);

% Plot median
v_median = median(all_variance, 1);
e_median = median(all_winderror, 1);

h_vmed = plot(h_vax, t_variance, v_median./numel(VW_grid), 'k-', 'linewidth', 1.5);
h_emed = plot(h_eax, t_variance, e_median, 'k-', 'linewidth', 1.5);

legend(h_vax, [h_var(2), h_var(1), h_vmed], {'Sample trajectory', 'Poorest mapping performance', 'Median'});
legend(h_eax, [h_err(2), h_err(1), h_emed], {'Sample trajectory', 'Poorest mapping performance', 'Median'});

