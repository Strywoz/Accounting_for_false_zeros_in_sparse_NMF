% For each sparsity S and each target density, run HALS variants on matrices of size m x m
% and plot the relative error vs lambda, for multiple random initializations.

% === Timestamp ===
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');

% === Parameters ===
lambda_list = [1, 0.5, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6];
densities_target_base = [0.9, 0.7, 0.5, 0.3, 0.1];
S_list = [0.9, 0.7, 0.5, 0.3, 0.1];
m = 200;

results = struct();

outdir = fullfile(pwd, sprintf('results_exp2_%s', timestamp));
if ~exist(outdir, 'dir'), mkdir(outdir); end

fprintf('\n=== Starting experiment 2 ===\n');

for si = 1:length(S_list)
    S = S_list(si);
    fprintf('\n=== Sparsity S=%.2f ===\n', S);

    density = 1 - S + 0.04;
    densities_target = densities_target_base(densities_target_base <= density);
    densities_target = sort(densities_target, 'descend');

    for di = 1:length(densities_target)
        density_target = densities_target(di);
        fprintf('\n--- Density Target %.0f%% ---\n', density_target*100);

        figure;
        hold on;
        colors = lines(5); % Generate distinct colors

        M2_all = zeros(5, length(lambda_list));

        for i = 1:5
            fprintf('\n=== initialisation %d/5 ===\n', i);

            final_m2 = run_hals_variants( ...
                m, m, 10, 1e10, 0.99, 300, ...
                false, false, true, S, [density_target], ...
                lambda_list, 42, randi(1e6));

            % Store ONLY relative error
            M2_all(i, :) = final_m2(1, :);

            % Plot
            idx = lambda_list > 0;
            loglog(lambda_list(idx), final_m2(1,idx), 's-', ...
                'LineWidth', 1.5, 'Color', colors(i,:));
        end

        % === Store all 5 initializations ===
        results.S(si).value = S;
        results.S(si).d(di).value = density_target;
        results.S(si).d(di).lambda = lambda_list;
        results.S(si).d(di).M2 = M2_all;


        xlabel('\lambda');
        ylabel('Relative Error (%)');
        title(sprintf('D=%.0f%% | S=%.0f%%', density_target*100, S*100));
        grid on;
        % legend('location', 'northwest');

        fname = sprintf('S%.1f_d%.1f.jpg', S, density_target);
        print(gcf, fullfile(outdir, fname), '-djpeg', '-r300');
        close(gcf);
    end
end

% === Final save ===
params = struct();
params.lambda_list = lambda_list;
params.S_list      = S_list;
params.m           = m;
params.timestamp   = timestamp;
params.densities_target_base = densities_target_base;

filename = sprintf('results_exp2_%s.mat', timestamp);
save(fullfile(outdir, filename), 'results', 'params');