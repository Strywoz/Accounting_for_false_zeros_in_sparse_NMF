% ============================================================
% Test HALS M2 on real data matrices
% ============================================================

clear; clc;

% === Parameters ===
lambda_list = [1, 0.5, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6];
densities_target = [0.9, 0.7, 0.5, 0.3, 0.1];

r = 10;
max_iter = 1e10;   
tol = 0.99;
max_time = 1e10;      

% === Input files ===
files = dir(fullfile('./data/', '*.mat'));

% === Output directory ===
outdir = fullfile(pwd, 'results_test_data');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

fprintf('\n=== Starting test on real data ===\n');

% ============================================================
% Loop over data files
% ============================================================
for i = 1:length(files)

    % ---- Load data ----
    file = fullfile(files(i).folder, files(i).name);
    data = load(file);

    if ~isfield(data, 'X')
        warning('File %s does not contain variable X — skipped.', files(i).name);
        continue;
    end

    X = data.X;

    % ---- Basic info ----
    real_density = nnz(X) / numel(X);
    fprintf('\n========================================\n');
    fprintf('File: %s\n', files(i).name);
    fprintf('Matrix size: %dx%d\n', size(X,1), size(X,2));
    fprintf('Real density: %.3f\n', real_density);
    fprintf('========================================\n');

    % ---- Run M2 ----
    [final_m2, W2_list, H2_list] = run_M2( ...
        X, r, max_iter, tol, max_time, ...
        lambda_list, densities_target);

    % ========================================================
    % Plot results
    % ========================================================
    figure;
    hold on;
    colors = lines(length(densities_target));

    for d = 1:length(densities_target)
        loglog(lambda_list, final_m2(d,:), 's-', ...
            'LineWidth', 1.5, ...
            'Color', colors(d,:), ...
            'DisplayName', sprintf('density = %.2f', densities_target(d)));
    end

    xlabel('\lambda');
    ylabel('Relative Error');
    title(sprintf('M2 on data: %s', files(i).name), 'Interpreter', 'none');
    grid on;
    legend('location', 'best');

    % ---- Save figure ----
    [~, base, ~] = fileparts(files(i).name);
    fig_name = sprintf('test_data_%s.jpg', base);
    print(gcf, fullfile(outdir, fig_name), '-djpeg', '-r300');
    close(gcf);

    % ========================================================
    % Save numerical results
    % ========================================================
    results = struct();
    results.lambda_list = lambda_list;
    results.densities_target = densities_target;
    results.final_m2 = final_m2;
    results.real_density = real_density;
    results.factors.W = W2_list;
    results.factors.H = H2_list;


    res_name = sprintf('results_test_%s.mat', base);
    save(fullfile(outdir, res_name), 'results');

end

fprintf('\n=== All tests completed ===\n');
