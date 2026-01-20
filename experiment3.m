% Percentage of true factorizations for M2 and NMF
% Compares the number of times M2 and NMF reach true factorization
% over multiple random matrices and initializations.

% === Timestamp ===
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');

% === Parameters ===
lambda_list = [1e-2];
densities_target_base = [0.9, 0.7, 0.5, 0.3, 0.1];
S_list = [0.9, 0.7, 0.5, 0.3, 0.1];

m = 200;

number_true_factorization_M2 = zeros(length(densities_target_base), length(S_list));
number_true_factorization_NMF = zeros(length(densities_target_base), length(S_list));
NMF_better_count = zeros(length(densities_target_base), length(S_list));

%% ============================================================
%  GLOBAL GRID OF DENSITIES
% ============================================================

densities_all = sort(densities_target_base, 'descend');
nD = length(densities_all);
nS = length(S_list);

outdir = fullfile(pwd, sprintf('results_exp3_%s', timestamp));
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% ============================================================
%  MAIN LOOP
% ============================================================

fprintf('\n=== Starting experiment 3 ===\n');

for si = 1:nS

    s = S_list(si);
    fprintf('\n=== Sparsity S = %.2f ===\n', s);

    density_max = 1 - s + 0.04;

    % Densities actually used for this S
    densities_target = densities_target_base(densities_target_base <= density_max);
    densities_target = sort(densities_target, 'descend');

    n_density = length(densities_target);

    % === HALS ===
    for n_matrix = 1:10

        fprintf('\n--- Matrix %d / 10 ---\n', n_matrix);
        seed = randi(1e6);
        for n_initialization = 1:10
            fprintf('Initialization %d / 10\n', n_initialization);
            [final_m0, final_m2] = run_hals_variants( ...
                m, m, 10, 1e10, 0.99, 1e1000, ...
                true, false, true, s, densities_target, lambda_list, seed, randi(1e6));
        
            for di = 1:n_density
                if final_m2(di,1) < 1e-1
                    number_true_factorization_M2(di, si) = number_true_factorization_M2(di, si) + 1;
                end
                if final_m0(di) < 1e-1
                    number_true_factorization_NMF(di, si) = number_true_factorization_NMF(di, si) + 1;
                end
                if final_m0(di) < final_m2(di, 1)
                    NMF_better_count(di, si) = NMF_better_count(di, si) + 1;
                end
            end
        end
    end
end

percent = struct();
percent.percent_true_factorization_M2 = number_true_factorization_M2
percent.percent_true_factorization_NMF = number_true_factorization_NMF


params = struct();
params.lambda_list = lambda_list;
params.S_list = S_list;
params.m = m;

% === final save ===

filename = sprintf('results_exp3_%s.mat', timestamp);
save(fullfile(outdir, filename), 'percent', 'NMF_better_count', 'params');