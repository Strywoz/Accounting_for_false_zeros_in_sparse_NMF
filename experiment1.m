% For each sparsity S and each target density, run HALS variants on matrices of different sizes
% and plot the relative error vs lambda.

% === Timestamp ===
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');

% === Parameters ===
lambda_list = [1, 0.5, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6];
densities_target_base = [0.9, 0.7, 0.5, 0.3, 0.1];
S_list = [0.9, 0.7, 0.5, 0.3, 0.1];
m_list = [100, 200, 500, 1000, 2000];

results = struct();

outdir = fullfile(pwd, sprintf('results_exp1_%s', timestamp));
if ~exist(outdir, 'dir'), mkdir(outdir); end

fprintf('\n=== Starting experiment 1 ===\n');

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
        colors = lines(length(m_list)); % Generate distinct colors

        for mi = 1:length(m_list)
            m = m_list(mi);
            fprintf('\n=== Matrix %dx%d ===\n', m, m);

            final_m2 = run_hals_variants( ...
                m, m, 10, 1e10, 0.99, 1e100, ...
                false, false, true, S, [density_target], lambda_list);

            % === Storage ===
            results.S(si).value = S;
            results.S(si).d(di).value = density_target;
            results.S(si).d(di).m(mi).m = m;
            results.S(si).d(di).m(mi).lambda = lambda_list;
            results.S(si).d(di).m(mi).M2 = final_m2;

            % === Plots ===
            idx = lambda_list > 0;
            loglog(lambda_list(idx), final_m2(1,idx), 's-', ...
                'LineWidth', 1.5, 'Color', colors(mi,:), ...
                'DisplayName', sprintf('m=%d', m));
        end

        xlabel('\lambda');
        ylabel('Relative Error (%)');
        title(sprintf('D=%.0f%% | S=%.0f%%', density_target*100, S*100));
        grid on;
        legend('location', 'southwest');

        fname = sprintf('S%.1f_d%.1f.jpg', S, density_target);
        print(gcf, fullfile(outdir, fname), '-djpeg', '-r300');
        close(gcf);
    end
end

% === Final save ===
params = struct();
params.lambda_list = lambda_list;
params.S_list      = S_list;
params.m_list      = m_list;
params.timestamp   = timestamp;
params.densities_target_base = densities_target_base;

filename = sprintf('results_exp1_%s.mat', timestamp);
save(fullfile(outdir, filename), 'results', 'params');