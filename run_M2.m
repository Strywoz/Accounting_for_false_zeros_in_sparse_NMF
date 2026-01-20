function varargout = run_M2(X_full, r, max_iter, tol, max_time, lambda_list, densities_target, rgn_seed2 = 43)
    % ===============================================================
    % Runs method 2.
    
    % Inputs :
    %   X_full     : full input matrix (m x n)
    %   r          : factorization rank
    %   max_iter   : maximum number of iterations
    %   tol        : stopping tolerance
    %   max_time   : maximum execution time (in seconds)
    %   lambda_list : list of lambda values to test
    %   densities_target  : list of densities to test
    %   rgn_seed2  : seed for generating the initializations W and H

    % Outputs :
    %   final_m2   : (n_density x n_lambda)
    %   W2_list   : cell array of W factors
    %   H2_list   : cell array of H factors
    % ===============================================================

    rng(rgn_seed2);
    rand_int = randi(1e6);

    % ======= true Matrix ========
    m = size(X_full, 1);
    n = size(X_full, 2);

    % ======= General Parameters ========
    % densities_target = [0.8, 0.6, 0.4, 0.2, 0.05];  % desired densities
    X_vec = sort(X_full(:), 'descend');
    N = m*n;

    density_list = zeros(size(densities_target));

    for i = 1:length(densities_target)
        k = max(floor(densities_target(i) * N),1);  % number of elements to keep
        density_list(i) = X_vec(k);                 % threshold for this density
    end
    
    n_density = length(density_list);
    n_lambda = length(lambda_list);

    % ======= Error storage ========
    final_m2 = zeros(n_density, n_lambda);

    W2_list = cell(n_density, n_lambda);
    H2_list = cell(n_density, n_lambda);

    start_time = tic;

    % =====================================================
    % Loop over densities
    % =====================================================
    for d = 1:n_density

        % rng('shuffle'); 
        rng(rgn_seed2); 
        rand_int = randi([1, 1000000]);

        threshold = density_list(d);

        X = X_full;
        X(X < threshold) = 0;

        fprintf("\n\nDensity %d : threshold = %.3f — actual density = %.3f\n\n", d, threshold, nnz(X)/(m*n));

        % Masks on true zeros and on values >0 for the observed matrix
        % M1 = (X_full == 0);
        % M2 = (X ~= 0);
        % M = M1 | M2;

        % Index of non-zero elements per row/column
        Zw_tilde = cell(m,1);
        Zh_tilde = cell(n,1);
        % Zw = cell(m,1);
        % Zh = cell(n,1);

        for i = 1:m
            Zw_tilde{i} = find(X(i,:) ~= 0);
            % Zw{i}       = find(X(i,:) == 0);
        end
        for j = 1:n
            Zh_tilde{j} = find(X(:,j) ~= 0);
            % Zh{j}       = find(X(:,j) == 0);
        end


        % =====================================================
        % LOOP OVER LAMBDAS (methods 1 and 2)
        % =====================================================
        for L = 1:n_lambda

            lambda = lambda_list(L);
            alpha = lambda - 1;

            % -----------------------------------------------------
            % METHOD 2 - quadratic version
            % -----------------------------------------------------
            fprintf("MMethod 2 (quadratic version) -- lambda = %.2e\n", lambda);

            rng(rand_int);  
            W2 = rand(m, r);
            H2 = rand(r, n);
            t = tic;
            prev_error = Inf;

            for iter = 1:max_iter

                % ---- W ----
                B = H2 * H2';
                A = X * H2';

                for q = 1:m
                    % col = Zw{q};
                    % B_tilde = H2(:, col)*H2(:, col)';
                    
                    col = Zw_tilde{q};
                    Hcol = H2(:, col);
                    B_tilde = B - Hcol * Hcol';

                    for p = 1:r
                        W2(q,p) = max(0, ...
                            (A(q,p) - W2(q,:)*B(:,p) + W2(q,p)*B(p,p) - alpha * (W2(q,:)*B_tilde(:,p) - W2(q,p)*B_tilde(p,p))) ...
                            / (B(p,p) + alpha*B_tilde(p,p) + 1e-10));
                    end
                end

                % ---- H ----
                B = W2' * W2;
                A = W2' * X;

                for j = 1:n
                    % row = Zh{j};
                    % B_tilde = W2(row,:)' * W2(row,:);

                    row = Zh_tilde{j};
                    W_row = W2(row, :);
                    B_tilde = B - W_row' * W_row;

                    for p = 1:r
                        H2(p,j) = max(0, ... 
                            (A(p,j) - B(p,:)*H2(:,j) + B(p,p)*H2(p,j) - alpha * (B_tilde(p,:)*H2(:,j) - B_tilde(p,p)*H2(p,j))) ...
                            / (B(p,p) + alpha*B_tilde(p,p) + 1e-10));
                    end
                end

                % ---- Critère ----
                R = X - W2*H2;
                current_error = norm(R,'fro')^2 + alpha * norm((W2*H2).*(X==0),'fro')^2;

                % error on full matrix
                final_m2(d,L) = norm(X_full - W2*H2,'fro') / norm(X_full,'fro');

                if toc(t) > max_time
                    fprintf("MMethod 2 λ=%.3f cut off at %.2f s (%.i iterations)\n", lambda, toc(t), iter);
                    break;
                end

                % if iter > 1 && abs(prev_error - current_error)/ (abs(prev_error) + 1) < tol
                %     fprintf("Cut off by tolerance tol = %.2e (%.i iterations)\n", tol, iter);
                %     break;
                % end

                if iter == 1
                    prev_error = current_error;
                end

                if mod(iter, 10) == 0
                    if prev_error * tol < current_error
                        fprintf("Cut off by tolerance tol = %.3f (%.i iterations in %.2f s)\n", tol, iter, toc(t));
                        break;
                    else
                        prev_error = current_error;
                    end
                end

                if iter == max_iter
                    fprintf("Reached the maximum number of iterations (%d) (in %.2f s)\n", max_iter, toc(t));
                end

                % prev_error = current_error;
            end
            W2_list{d,L} = W2;
            H2_list{d,L} = H2;
        end
    end

    fprintf("\nTotal time: %.2f s\n", toc(start_time));

    % ===========================================================
    % Prepare output according to booleans
    % ===========================================================

    % afficher les erreurs finales en %
    out = {}; 
    out{end+1} = final_m2;
    out{end+1} = W2_list;
    out{end+1} = H2_list;

    varargout = out;
end