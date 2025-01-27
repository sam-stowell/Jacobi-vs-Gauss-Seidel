function Assignment_graphs(n, abc)
    clc;

    % Matrix sizes to test
    matrix_sizes = [10, 50, 100, 200, 500, 1000, 5000, 10000, 20000];


    % Preallocate arrays to store times
    jacobi_times = zeros(length(matrix_sizes), 1);
    gauss_seidel_times = zeros(length(matrix_sizes), 1);
    jacobi_setup_times = zeros(length(matrix_sizes), 1);
    gauss_seidel_setup_times = zeros(length(matrix_sizes), 1);
    jacobi_total_times = zeros(length(matrix_sizes), 1);
    gauss_seidel_total_times = zeros(length(matrix_sizes), 1);

    % Loop through each matrix size
    for i = 1:length(matrix_sizes)
        n = matrix_sizes(i);

        % Measure total time for Jacobi method
        disp('Running Jacobi');
        disp(n);
        [jacobi_computation_time, jacobi_setup_time, iter_jacobi] = jacobi_method(n, abc);
        jacobi_times(i) = jacobi_computation_time;
        jacobi_setup_times(i) = jacobi_setup_time;
        jacobi_total_times(i) = jacobi_computation_time + jacobi_setup_time;

        % Measure total time for Gauss-Seidel method
        disp('Running Gauss-Seidel');
        disp(n);
        [gauss_seidel_computation_time, gauss_seidel_setup_time, iter_gs] = gauss_seidel_method(n, abc);
        gauss_seidel_times(i) = gauss_seidel_computation_time;
        gauss_seidel_setup_times(i) = gauss_seidel_setup_time;
        gauss_seidel_total_times(i) = gauss_seidel_computation_time + gauss_seidel_setup_time;
    end

    % Create figure with 3 subplots
    figure;

    % Subplot 1 Computation Time
    % 3 rows, 1 column first subplot
    subplot(3, 1, 1);
    % Plot computation time vs matrix sizes jacobi
    plot(matrix_sizes, jacobi_times, '-o', 'LineWidth', 1.5, 'DisplayName', 'Jacobi Method');
    % Plot added to current figure
    hold on;
    % Plot computation time vs matrix sizes GS
    plot(matrix_sizes, gauss_seidel_times, '-x', 'LineWidth', 1.5, 'DisplayName', 'Gauss-Seidel Method');
    % X axis label
    xlabel('Matrix Size (n)');
    % Y axis label
    ylabel('Computation Time (seconds)');
    % Title
    title('Computation Time vs Matrix Size');
    legend('Location', 'northwest');
    grid on;

    % Subplot 2 Matrix Setup Time
    % 3 rows, 1 column second subplot
    subplot(3, 1, 2);
    % Plot setup time vs matrix sizes jacobi
    plot(matrix_sizes, jacobi_setup_times, '-o', 'LineWidth', 1.5, 'DisplayName', 'Jacobi Method');
    % Plot added to current figure
    hold on;
    % Plot setup time vs matrix sizes GS
    plot(matrix_sizes, gauss_seidel_setup_times, '-x', 'LineWidth', 1.5, 'DisplayName', 'Gauss-Seidel Method');
    % X axis label
    xlabel('Matrix Size (n)');
    % Y axis label
    ylabel('Matrix Setup Time (seconds)');
    % Title
    title('Matrix Setup Time vs Matrix Size');
    legend('Location', 'northwest');
    grid on;

    % Subplot 3 Total Time
    % 3 rows, 1 column third subplot
    subplot(3, 1, 3);  
    % Plot total time vs matrix sizes jacobi
    plot(matrix_sizes, jacobi_total_times, '-o', 'LineWidth', 1.5, 'DisplayName', 'Jacobi Method');
    % Plot added to current figure
    hold on;
    % Plot total time vs matrix sizes jacobi
    plot(matrix_sizes, gauss_seidel_total_times, '-x', 'LineWidth', 1.5, 'DisplayName', 'Gauss-Seidel Method');
    % X axis label
    xlabel('Matrix Size (n)');
    % Y axis label
    ylabel('Total Time (seconds)');
    % Title
    title('Total Time vs Matrix Size');
    legend('Location', 'northwest');
    grid on;

    % Display results in tabular form
    results_table = table(matrix_sizes', jacobi_setup_times, jacobi_times, jacobi_total_times, gauss_seidel_setup_times, gauss_seidel_times, gauss_seidel_total_times,'VariableNames', {'Matrix_Size', 'Jacobi_Setup_Time', 'Jacobi_Computation_Time', 'Jacobi_Total_Time', 'Gauss_Seidel_Setup_Time', 'Gauss_Seidel_Computation_Time', 'Gauss_Seidel_Total_Time'});

    % Display the table
    disp('Computation Time for Jacobi vs Gauss-Seidel:');
    disp(results_table);
end
% 
%     clc;
% 
%     % Matrix size limit
%     matrix_max = 1000;  % Change this to any value up to which you want to test
% 
%     % Preallocate arrays to store times
%     jacobi_times = zeros(matrix_max, 1);  % Preallocate for all sizes up to matrix_max
%     gauss_seidel_times = zeros(matrix_max, 1);
%     jacobi_setup_times = zeros(matrix_max, 1);
%     gauss_seidel_setup_times = zeros(matrix_max, 1);
%     jacobi_total_times = zeros(matrix_max, 1);
%     gauss_seidel_total_times = zeros(matrix_max, 1);
% 
%     % Loop through each matrix size from 1 to matrix_max
%     for n = 1:matrix_max
%         % Measure total time for Jacobi method
%         disp('Running Jacobi');
%         disp(n);
%         [jacobi_computation_time, jacobi_setup_time, iter_jacobi] = jacobi_method(n, abc);
%         jacobi_times(n) = jacobi_computation_time;
%         jacobi_setup_times(n) = jacobi_setup_time;
%         jacobi_total_times(n) = jacobi_computation_time + jacobi_setup_time;
% 
%         % Measure total time for Gauss-Seidel method
%         disp('Running Gauss-Seidel');
%         disp(n);
%         [gauss_seidel_computation_time, gauss_seidel_setup_time, iter_gs] = gauss_seidel_method(n, abc);
%         gauss_seidel_times(n) = gauss_seidel_computation_time;
%         gauss_seidel_setup_times(n) = gauss_seidel_setup_time;
%         gauss_seidel_total_times(n) = gauss_seidel_computation_time + gauss_seidel_setup_time;
%     end
% 
%     % Create figure with 3 subplots
%     figure;
% 
%     % Subplot 1 Computation Time
%     % 3 rows, 1 column first subplot
%     subplot(3, 1, 1);
%     % Plot computation time vs matrix sizes jacobi
%     plot(1:matrix_max, jacobi_times, '-o', 'LineWidth', 1.5, 'DisplayName', 'Jacobi Method');
%     % Plot added to current figure
%     hold on;
%     % Plot computation time vs matrix sizes GS
%     plot(1:matrix_max, gauss_seidel_times, '-x', 'LineWidth', 1.5, 'DisplayName', 'Gauss-Seidel Method');
%     % X axis label
%     xlabel('Matrix Size (n)');
%     % Y axis label
%     ylabel('Computation Time (seconds)');
%     % Title
%     title('Computation Time vs Matrix Size');
%     legend('Location', 'northwest');
%     grid on;
% 
%     % Subplot 2 Matrix Setup Time
%     % 3 rows, 1 column second subplot
%     subplot(3, 1, 2);
%     % Plot setup time vs matrix sizes jacobi
%     plot(1:matrix_max, jacobi_setup_times, '-o', 'LineWidth', 1.5, 'DisplayName', 'Jacobi Method');
%     % Plot added to current figure
%     hold on;
%     % Plot setup time vs matrix sizes GS
%     plot(1:matrix_max, gauss_seidel_setup_times, '-x', 'LineWidth', 1.5, 'DisplayName', 'Gauss-Seidel Method');
%     % X axis label
%     xlabel('Matrix Size (n)');
%     % Y axis label
%     ylabel('Matrix Setup Time (seconds)');
%     % Title
%     title('Matrix Setup Time vs Matrix Size');
%     legend('Location', 'northwest');
%     grid on;
% 
%     % Subplot 3 Total Time
%     % 3 rows, 1 column third subplot
%     subplot(3, 1, 3);  
%     % Plot total time vs matrix sizes jacobi
%     plot(1:matrix_max, jacobi_total_times, '-o', 'LineWidth', 1.5, 'DisplayName', 'Jacobi Method');
%     % Plot added to current figure
%     hold on;
%     % Plot total time vs matrix sizes jacobi
%     plot(1:matrix_max, gauss_seidel_total_times, '-x', 'LineWidth', 1.5, 'DisplayName', 'Gauss-Seidel Method');
%     % X axis label
%     xlabel('Matrix Size (n)');
%     % Y axis label
%     ylabel('Total Time (seconds)');
%     % Title
%     title('Total Time vs Matrix Size');
%     legend('Location', 'northwest');
%     grid on;
% 
%     % Display results in tabular form
%     results_table = table((1:matrix_max)', jacobi_setup_times, jacobi_times, jacobi_total_times, ...
%                           gauss_seidel_setup_times, gauss_seidel_times, gauss_seidel_total_times, ...
%                           'VariableNames', {'Matrix_Size', 'Jacobi_Setup_Time', 'Jacobi_Computation_Time', ...
%                                             'Jacobi_Total_Time', 'Gauss_Seidel_Setup_Time', 'Gauss_Seidel_Computation_Time', 'Gauss_Seidel_Total_Time'});
% 
%     % Display the table
%     disp('Computation Time for Jacobi vs Gauss-Seidel:');
%     disp(results_table);
% end


% Function generate matrix with the output being A
function A = generate_matrix(n)
    % Create a random matrix and make it strictly diagonally dominant
    A = randi([-10,10],n);
    for i = 1:n
        % Ensure strict diagonal dominance
        % A(i, i) = sum(abs(A(i, :))) + randi([-10,10]);
        A(i, i) = sum(abs(A(i, :))) - abs(A(i, i)) + randi([1, 10]);
    end
end

% To return computation_time, setup_time, iter for the graphs
function [computation_time, setup_time, iter] = jacobi_method(n, abc)
    % Randomly generated data with abc as last 3 non 0 digits in student number
    % 447
    rng(abc);

    % Start time for matrix setup
    tic;
    % Generate matrix A
    A = generate_matrix(n);
    % Generate random right-hand side vector b
    b = randi([-10, 10], n, 1);

    % Convergence tolerance stated in assignment brief increased instead of 0.001
    tol = 0.00001;
    % Maximum iterations
    max_iter = 1000;

    % Ax = b
    % A - Coefficient matrix
    % b - Right-hand side vector
    % x0 - Initial guess vector
    % tol - Convergence tolerance
    % max_iter - Maximum number of iterations
    % x - Solution vector
    % iter - Number of iterations taken

    
    % Diagonal matrix (diagonal elements of A)
    D = diag(diag(A));
    % Lower
    L = tril(A, -1);
    % Upper
    U = triu(A, 1);

    % Initial guess
    x0 = ones(n, 1); 
    % End time for matrix setup
    setup_time = toc; 
    
    disp('Run Jacobi Method')
    % Start timing computation
    tic; 
    % Run Jacobi Method
    [x, iter,] = jacobi_iteration(L, D, U, x0, b, tol, max_iter);
    % End timing computation
    computation_time = toc;

    % Display results
    disp('Results:')
    disp('x:')
    disp(x);
    disp('Matrix Setup Time (seconds): ')
    disp(setup_time)
    disp('Jacobi Method: Iterations')
    disp(iter)
    disp('Time (seconds):')
    disp(computation_time)
    disp('Total time:')
    disp(computation_time+setup_time)
end

% Function generate matrix with the output being x and iterations
function [x, iter] = jacobi_iteration(L, D, U, x0, b, tol, max_iter)
    
    x=x0;
    % Jacobi iteration
    for iter = 1:max_iter
        % Jacobi update formula including inverse D
        x_new = D \ (b - (L + U) * x);

        % Check for convergence
        % If the maximum absolute difference is smaller than tolerance
        if norm(x_new - x, inf) < tol
            x = x_new;
            return;
        end

        % Update the solution
        x = x_new;
    end
    
    % Maximum iterations reached without convergence
    error('Jacobi method did not converge within the maximum number of iterations.');
end

% To return computation_time, setup_time, iter for the graphs
function [computation_time, setup_time, iter] = gauss_seidel_method(n, abc)
    % Randomly generated data with abc as last 3 non 0 digits in student number
    % 447
    rng(abc);

    % Start time
    tic; 
    % Generate matrix A
    A = generate_matrix(n);
    % Generate random right-hand side vector b
    b = randi([-10, 10], n, 1);

    % Convergence tolerance stated in the assignment breif increased instead of 0.001
    tol = 0.00001;
    % Maximum number of iterations
    max_iter = 1000;

    % Ax = b
    % A - Coefficient matrix (n x n)
    % b - Right-hand side vector (n x 1)
    % x0 - Initial guess vector (n x 1)
    % tol - Convergence tolerance
    % max_iter - Maximum number of iterations
    % x - Solution vector (n x 1)
    % iter - Number of iterations performed

    % Diagonal matrix (diagonal elements of A)
    D = diag(diag(A));
    % Lower
    L = tril(A, -1);
    % Upper
    U = triu(A, 1);
    
    % Initial guess
    n = length(b);
    x0 = ones(n, 1); 
    % End time for matrix setup
    setup_time = toc; 

    disp('Run Jacobi Method')
    % Start timing computation
    tic; 
    % Run GS Method
    [x, iter,] = gauss_seidel_iteration(L, D, U, x0, b, tol, max_iter);
    % End timing computation
    computation_time = toc;

    % Display results
    disp('Results:')
    disp('x:')
    disp(x);
    disp('Matrix Setup Time (seconds): ')
    disp(setup_time)
    disp('Gauss Seidel Method: Iterations')
    disp(iter)
    disp('Time (seconds):')
    disp(computation_time)
    disp('Total time:')
    disp(computation_time+setup_time)
end

function [x, iter] = gauss_seidel_iteration(L, D, U, b, x0, tol, max_iter)

    x=x0;
    % Gauss-Seidel iteration
    for iter = 1:max_iter
        
        x_new = (L+D) \ (b-(U*x0));
        if norm((x_new - x), inf) < tol
            x = x_new;
            return;
        end
        % Update the solution
        x = x_new;
    end
    
    % Maximum iterations reached without convergence
    error('GS method did not converge within the maximum number of iterations.');
end