% Plot result
plot_proportion = 250/500;

recalculate_GP = false;
if ~exist('GPt', 'var'); GPt = false; end


figure(6); clf;
set(gca, 'Zdir', 'reverse'); set(gca, 'Ydir', 'reverse'); axis equal;
view(3); hold on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
h_box = plot3(x_limits([1,2,2,1,1,1,2,2,1,1,1,1,2,2,2,2]), ...
    y_limits([1,1,1,1,1,2,2,2,2,2,2,1,1,2,2,1]), ...
    z_limits([1,1,2,2,1,1,1,2,2,1,2,2,2,2,1,1]), ...
    'Color', [.8, .8, .8], 'LineStyle', '--');

h_start = plot3(start_pos(1), start_pos(2), start_pos(3), 'g^', 'Color',...
    [0,.5,0],'MarkerSize', 10, 'LineWidth', 2);

h_preflight = plot3(X_init(:,1), X_init(:,2), X_init(:,3), 'g-');

% h_energy_target = plot3(energy_target(1), energy_target(2), ...
%     energy_target(3), 'x','Color', [0.5,0.5,1], 'MarkerSize', 10);
% 
% h_info_target = plot3(info_target(1), info_target(2), ...
%     info_target(3), 'bx','Color', [1,0.5,0], 'MarkerSize', 10);

CTrain = squeeze(X_train_full(:,:,floor(n_replan*plot_proportion)));
WTrain = squeeze(W_train_full(:,:,floor(n_replan*plot_proportion)));
if GPt; tTrain = t_train_full(:,floor(n_replan*plot_proportion)); end;
cc = max_obs; while CTrain(cc,1)==0; cc=cc-1; end;
CTrain = CTrain(1:cc, :); WTrain = WTrain(1:cc, :); 
h_train = plot3(CTrain(:,1), CTrain(:,2), CTrain(:,3), 'r+', 'MarkerSize', 3);

last_index = floor(size(pos_full, 2)*plot_proportion);
h_path = plot3(pos_full(1,1:last_index), pos_full(2,1:last_index), ...
    pos_full(3,1:last_index), 'b-');

if exist('thermals', 'var')
	plot3(thermals(:,4), thermals(:,5), thermals(:,6), 'o', ...
		'color', [1, 0.6, 0.1], 'markersize', 10, 'LineWidth', 2, ...
		'MarkerFaceColor', [1, 0.6, 0.1]);
end
		
	

axis tight;

%% Recalculate GP

if recalculate_GP
opt = optimset('Display', 'iter', 'GradObj', 'on', 'TolX', 1e-7,...
    'MaxIter', 100, 'maxfunevals', 2000);

if GPt
	[loghyper2, fmin] = fminunc(@(hyper) ...
		GPt_likelihoodn(CTrain, tTrain, WTrain(:,1), WTrain(:,2), WTrain(:,3), cov_funs, hyper), loghyper, opt);
	fprintf(1, '\nRetrain: l = %0.5g, lt = %0.5g, sf = %0.5g, sn = %0.5g\n', exp(loghyper2));
% 	fprintf(1, ['\nRetrain: l = %0.5g, lt = %0.5g, ld = %0.5g sf = %0.5g, \n', ...
% 		'K = %0.5g, Wx = %0.5g, Wy = %0.5g, Wz = %0.5g, sn = %0.5g\n'], ...
% 		[exp(loghyper2(1:4)), loghyper2(5:8), exp(loghyper2(9))]);
	
else	
	[loghyper2, fmin] = fminunc(@(hyper) ...
		GP_likelihoodn(CTrain, WTrain(:,1), WTrain(:,2), WTrain(:,3), cov_funs, hyper), loghyper, opt);
	fprintf(1, '\nRetrain: l = %0.5g, sf = %0.5g, sn = %0.5g\n', exp(loghyper2));
end

elseif exist('save_hypers', 'var')
	save_row = sum(save_hypers(:,1) < floor(n_replan*plot_proportion));
	loghyper2 = save_hypers(save_row, 2:end);
end

%% Field estimate

figure(7); clf;
set(gca, 'Zdir', 'reverse'); set(gca, 'Ydir', 'reverse'); axis equal;
view(3); hold on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
h_box = plot3(x_limits([1,2,2,1,1,1,2,2,1,1,1,1,2,2,2,2]), ...
    y_limits([1,1,1,1,1,2,2,2,2,2,2,1,1,2,2,1]), ...
    z_limits([1,1,2,2,1,1,1,2,2,1,2,2,2,2,1,1]), ...
    'Color', [.8, .8, .8], 'LineStyle', '--');

h_train2 = plot3(CTrain(:,1), CTrain(:,2), CTrain(:,3), 'r+', 'MarkerSize', 3);

if GPt
	tnow = (plot_proportion*(tf-t0)+t0)*ones(numel(x_grid), 1);
	[WW_grid, VW_grid] = GPt_predict(CTrain, tTrain, WTrain, ...
		[x_grid(:), y_grid(:), z_grid(:)], tnow, cov_funs{1}, loghyper2);

else
	[WW_grid, VW_grid] = GP_predict(CTrain, WTrain, ...
		[x_grid(:), y_grid(:), z_grid(:)], cov_funs{1}, loghyper2);
end

U_grid = reshape(WW_grid(:,1), size(x_grid));
V_grid = reshape(WW_grid(:,2), size(x_grid));
W_grid = reshape(WW_grid(:,3), size(x_grid));
VW_grid = reshape(VW_grid, size(x_grid));
scale = sqrt(max(sum(WW_grid.^2, 2)))/max_wind*0.6;
h_windcone = coneplot(x_grid, y_grid, z_grid, U_grid, V_grid, W_grid, ...
    x_grid, y_grid, z_grid, VW_grid, scale);
set(h_windcone, 'EdgeColor', 'none')
axis tight
% h_cbar = colorbar('North')


%% Magicness

axlim = zeros(3,2);
axlim(1,:) = get(gca, 'XLim');
axlim(2,:) = get(gca, 'YLim');
axlim(3,:) = get(gca, 'ZLim');

%% Actual wind field
figure(8); clf;
t_now = plot_proportion*(tf-t0)+t0;
set(gca, 'Zdir', 'reverse'); set(gca, 'Ydir', 'reverse'); axis equal;
view(3); hold on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
h_box = plot3(x_limits([1,2,2,1,1,1,2,2,1,1,1,1,2,2,2,2]), ...
    y_limits([1,1,1,1,1,2,2,2,2,2,2,1,1,2,2,1]), ...
    z_limits([1,1,2,2,1,1,1,2,2,1,2,2,2,2,1,1]), ...
    'Color', [.8, .8, .8], 'LineStyle', '--');

[xx_grid, yy_grid, zz_grid] = meshgrid( ...
	linspace(x_limits(1),  x_limits(2), 8), ...
	linspace(y_limits(1),  y_limits(2), 8), ...
	linspace(z_limits(1),  z_limits(2), 5));


% if GPt
	[WA_grid] = W_actual([xx_grid(:), yy_grid(:), zz_grid(:)]', ...
		t_now*ones(1,numel(x_grid)))';
% else
% 	[WA_grid] = W_actual([xx_grid(:), yy_grid(:), zz_grid(:)]')';
% end
UU_grid = reshape(WA_grid(:,1), size(xx_grid));
VV_grid = reshape(WA_grid(:,2), size(xx_grid));
WW_grid = reshape(WA_grid(:,3), size(xx_grid));
h_windcone2 = coneplot(xx_grid, yy_grid, zz_grid, UU_grid, VV_grid, WW_grid, ...
    xx_grid, yy_grid, zz_grid, sqrt(max(sum(WA_grid.^2, 2)))/max_wind*0.6);
set(h_windcone2, 'EdgeColor', 'none', 'FaceColor', 'b')
axis tight

%% 
set(gca, 'XLim', axlim(1,:));
set(gca, 'YLim', axlim(2,:));
set(gca, 'ZLim', axlim(3,:));