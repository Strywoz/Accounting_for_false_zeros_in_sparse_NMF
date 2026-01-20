function varargout = run_hals_variants(m, n, r, max_iter, tol, max_time, doNMF, doM1, doM2, S, densities_target, lambda_list, rgn_seed1 = 42, rgn_seed2 = 43)
    % ===============================================================
    % Runs classic HALS method, method 1 and/or method 2.
    
    % Inputs :
    %   m, n, r    : NMF dimensions
    %   max_iter   : maximum number of iterations
    %   tol        : stopping tolerance
    %   max_time   : maximum execution time (in seconds)
    %   doNMF      : run method 0 (classic HALS)
    %   doM1       : run method 1 (linear)
    %   doM2       : run method 2 (quadratic)
    %   S          : target density for the full matrix (between 0 and 1)
    %   densities_target : vector of target densities (between 0 and 1)
    %   lambda_list : list of lambda values to test
    %   rgn_seed1  : seed for generating the full matrix
    %   rgn_seed2  : seed for generating the initializations W and H

    % Outputs (varargout) :
    %   final_m0   : (n_density x 1)
    %   final_m1   : (n_density x n_lambda)
    %   final_m2   : (n_density x n_lambda)
    % ===============================================================

    rng(rgn_seed1);

    % ======= true Matrix ========
    W_true = rand(m, r);
    H_true = rand(r, n);

    p = 1 - exp(log(1-exp(log(S)/10))/2);

    W_true(W_true < p) = 0;
    H_true(H_true < p) = 0;

    X_full = W_true * H_true;
    fprintf("\n\nDensity %d \n\n", nnz(X_full)/(m*n));

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
    final_m0 = zeros(n_density, 1);
    final_m1 = zeros(n_density, n_lambda);
    final_m2 = zeros(n_density, n_lambda);

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
        M1 = (X_full == 0);
        M2 = (X ~= 0);
        M = M1 | M2;

        % Index of non-zero elements per row/column
        Zw_tilde = cell(m,1);
        Zh_tilde = cell(n,1);
        Zw = cell(m,1);
        Zh = cell(n,1);

        for i = 1:m
            Zw_tilde{i} = find(X(i,:) ~= 0);
            Zw{i}       = find(X(i,:) == 0);
        end
        for j = 1:n
            Zh_tilde{j} = find(X(:,j) ~= 0);
            Zh{j}       = find(X(:,j) == 0);
        end

        % =====================================================
        % METHOD 0 : classic HALS
        % =====================================================
        if doNMF
            fprintf("MMethod 0 (classic HALS)\n");

            rng(rand_int);
            W0 = rand(m, r);
            H0 = rand(r, n);

            % Scaling
            A = X * H0';
            B = H0 * H0';
            scale = sum(sum(A.*W0)) / sum(sum(B.*(W0'*W0)));
            H0 = H0 * scale;

            t = tic;
            prev_error = Inf;

            for iter = 1:max_iter
                
                % ---- Update W ----
                B = H0 * H0';
                A = X * H0';

                for p = 1:r
                    W0(:,p) = max(0, (A(:,p) - W0*B(:,p) + W0(:,p)*B(p,p)) / (B(p,p) + 1e-10));
                end

                % ---- Update H ----
                B = W0' * W0;
                A = W0' * X;

                for p = 1:r
                    H0(p,:) = max(0, (A(p,:) - B(p,:)*H0 + B(p,p)*H0(p,:)) / (B(p,p) + 1e-10));
                end

                % ---- Criterion ----
                R = X - W0*H0;
                current_error = norm(R,'fro')^2;

                final_m0(d) = norm(R(M),'fro') / norm(X,'fro');

                if toc(t) > max_time
                    fprintf("MMethod 0 cut off at %.2f s (%.i iterations)\n", toc(t), iter);
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
        end

        % =====================================================
        % LOOP OVER LAMBDAS (methods 1 and 2)
        % =====================================================
        for L = 1:n_lambda

            lambda = lambda_list(L);
            alpha = lambda - 1;

            % -----------------------------------------------------
            % METHOD 1 - linear version
            % -----------------------------------------------------
            if doM1
                fprintf("MMethod 1 (linear version) -- lambda = %.2e\n", lambda);

                rng(rand_int);
                W1 = rand(m, r);
                H1 = rand(r, n);

                % Scaling
                A = X * H1';
                B = H1 * H1';
                scale = sum(sum(A.*W1)) / sum(sum(B.*(W1'*W1)));
                H1 = H1 * scale;

                t = tic;
                prev_error = Inf;

                for iter = 1:max_iter

                    % ---- W ----
                    B = H1 * H1';
                    A = X * H1';
                    S = sum(H1, 2);

                    for p = 1:r
                        som_H_col = zeros(m,1);
                        for i = 1:m
                            som_H_col(i) = S(p) - sum(H1(p, Zw_tilde{i}));
                        end
                        W1(:,p) = max(0, (A(:,p) - W1*B(:,p) + W1(:,p)*B(p,p) - (alpha/2)*som_H_col) / (B(p,p) + 1e-10));
                    end

                    % ---- H ----
                    B = W1' * W1;
                    A = W1' * X;
                    S = sum(W1,1);

                    for p = 1:r
                        som_W_row = zeros(1,n);
                        for j = 1:n
                            som_W_row(j) = S(p) - sum(W1(Zh_tilde{j}, p));
                        end
                        H1(p,:) = max(0, (A(p,:) - B(p,:)*H1 + B(p,p)*H1(p,:) - (alpha/2)*som_W_row) / (B(p,p) + 1e-10));
                    end

                    % ---- Criterion ----
                    R = X - W1*H1;
                    current_error = norm(R, 'fro')^2 + alpha * sum(sum((W1*H1).*(X==0)));

                    final_m1(d,L) = norm(R(M),'fro') / norm(X,'fro');

                    if toc(t) > max_time
                        fprintf("MMethod 1 λ=%.3f cut off at %.2f s (%.i iterations)\n", lambda, toc(t), iter);
                        break;
                    end

                    % stopping criterion different because "error" on the linear method can be negative !!
                    if iter > 1 && abs(prev_error - current_error)/ (abs(prev_error) + 1) < 1e-6
                        fprintf("Cut off by tolerance tol = %.3f (%.i iterations in %.2f s)\n", 1e-6, iter, toc(t));
                        break;
                    end

                    if iter == max_iter
                        fprintf("Reached the maximum number of iterations (%d) (in %.2f s)\n", max_iter, toc(t));
                    end

                    prev_error = current_error;
                end
            end

            % -----------------------------------------------------
            % METHOD 2 - quadratic version
            % -----------------------------------------------------
            if doM2
                fprintf("MMethod 2 (quadratic version) -- lambda = %.2e\n", lambda);

                rng(rand_int);  
                W2 = rand(m, r);
                H2 = rand(r, n);

                % Scaling
                A = X * H2';
                B = H2 * H2';
                scale = sum(sum(A.*W2)) / sum(sum(B.*(W2'*W2)));
                H2 = H2 * scale;

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

                    final_m2(d,L) = norm(R(M),'fro') / norm(X,'fro');

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
            end
        end
    end

    fprintf("\nTotal time: %.2f s\n", toc(start_time));

    % ===========================================================
    % Prepare output according to booleans
    % ===========================================================

    % afficher les erreurs finales en %
    out = {}; 
    if doNMF, out{end+1} = final_m0 * 100; end
    if doM1, out{end+1} = final_m1 * 100; end
    if doM2, out{end+1} = final_m2 * 100; end

    varargout = out;

end