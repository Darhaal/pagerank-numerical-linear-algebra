%% Application of Numerical Linear Algebra to PageRank
% MTH 4311 - Numerical Analysis
% Artem Okhten

clear;
clc;
close all;

%% Setup

n = 10;
alpha = 0.85;
TOL = 1e-10;
maxIter = 1000;

v = ones(n, 1) / n;

fprintf("n = %d, alpha = %.2f, tolerance = %.1e\n\n", n, alpha, TOL);

%% Sparse Web Graph

% A(i,j)=1 means page j links to page i.
% Page 10 is left as a dangling node.

row_idx = [ ...
    2, 3, 4, ...
    3, 5, ...
    1, 6, 7, ...
    3, 7, 8, ...
    4, 8, ...
    3, 9, ...
    1, 6, 10, ...
    5, 7, 10, ...
    6, 10 ...
];

col_idx = [ ...
    1, 1, 1, ...
    2, 2, ...
    3, 3, 3, ...
    4, 4, 4, ...
    5, 5, ...
    6, 6, ...
    7, 7, 7, ...
    8, 8, 8, ...
    9, 9 ...
];

A = sparse(row_idx, col_idx, 1, n, n);

fprintf("Nonzero links in A: %d\n", nnz(A));

%% Transition Matrix

[S, dangling, outDegree] = build_sparse_transition_matrix(A);

columnSums = full(sum(S, 1));
nonDanglingColumns = ~dangling;
maxColumnError = max(abs(columnSums(nonDanglingColumns) - 1));

fprintf("Nonzero entries in S: %d\n", nnz(S));
fprintf("Dangling pages: ");
disp(find(dangling));
fprintf("Max column-sum error: %.3e\n\n", maxColumnError);

%% Graph Figure

G = digraph(A');

figure;
p = plot(G, ...
    'Layout', 'circle', ...
    'NodeLabel', compose('P%d', 1:n));

title('Directed Web Graph');
p.MarkerSize = 8;
p.LineWidth = 1.2;

saveas(gcf, 'pagerank_graph.png');

%% Power Iteration

[r_power, errors, massErrors, numIter] = pagerank_power_sparse( ...
    S, dangling, alpha, v, TOL, maxIter);

fprintf("Power iteration iterations: %d\n", numIter);
fprintf("Final L1 change: %.3e\n", errors(end));
fprintf("Sum of PageRank vector: %.12f\n", sum(r_power));
fprintf("Max mass error before renormalization: %.3e\n\n", max(massErrors));

%% Results Table

pageNumbers = (1:n)';
inDegree = full(sum(A, 2));
outDegreeColumn = outDegree(:);

resultsTable = table( ...
    pageNumbers, ...
    inDegree, ...
    outDegreeColumn, ...
    dangling(:), ...
    r_power, ...
    'VariableNames', {'Page', 'InDegree', 'OutDegree', 'IsDangling', 'PageRank'} ...
);

rankedTable = sortrows(resultsTable, 'PageRank', 'descend');

disp("PageRank results:");
disp(rankedTable);

%% PageRank Values Figure

figure;
bar(1:n, r_power);
xlabel('Page Number');
ylabel('PageRank Value');
title(sprintf('Computed PageRank Values, alpha = %.2f', alpha));
grid on;

saveas(gcf, 'pagerank_values.png');

%% Convergence Figure

figure;
semilogy(1:length(errors), errors, '-o', 'LineWidth', 1.2);
xlabel('Iteration');
ylabel('||r^{(k+1)} - r^{(k)}||_1');
title(sprintf('Power Iteration Convergence, alpha = %.2f', alpha));
grid on;

saveas(gcf, 'pagerank_convergence.png');

%% Convergence Slope Estimate

iterIndex = (1:length(errors))';

% Skip first few points because early iterations are less linear.
fitStart = 5;

p = polyfit(iterIndex(fitStart:end), log10(errors(fitStart:end)), 1);

convSlope = p(1);
convFactor = 10^convSlope;

fprintf("Convergence slope on log10 scale: %.4f\n", convSlope);
fprintf("Approximate error factor per iteration: %.4f\n\n", convFactor);

%% Damping Factor Comparison

alphaValues = [0.70, 0.85, 0.95];

allErrors = cell(length(alphaValues), 1);
allIterations = zeros(length(alphaValues), 1);
finalErrors = zeros(length(alphaValues), 1);

for i = 1:length(alphaValues)
    currentAlpha = alphaValues(i);

    [~, errHist, ~, iters] = pagerank_power_sparse( ...
        S, dangling, currentAlpha, v, TOL, maxIter);

    allErrors{i} = errHist;
    allIterations(i) = iters;
    finalErrors(i) = errHist(end);
end

alphaTable = table(alphaValues(:), allIterations(:), finalErrors(:), ...
    'VariableNames', {'Alpha', 'Iterations', 'FinalError'});

disp("Damping factor comparison:");
disp(alphaTable);

figure;
hold on;

for i = 1:length(alphaValues)
    semilogy(1:length(allErrors{i}), allErrors{i}, ...
        'LineWidth', 1.4, ...
        'DisplayName', sprintf('\\alpha = %.2f', alphaValues(i)));
end

xlabel('Iteration');
ylabel('||r^{(k+1)} - r^{(k)}||_1');
title('Effect of Damping Factor on Convergence');
legend('Location', 'northeast');
grid on;
hold off;

saveas(gcf, 'pagerank_alpha_comparison.png');

%% Direct Solve Check

S_full = full(S);

for j = 1:n
    if dangling(j)
        S_full(:, j) = v;
    end
end

I = eye(n);
r_direct = (I - alpha*S_full) \ ((1 - alpha) * v);

directMassError = abs(sum(r_direct) - 1);
directDifference = norm(r_direct - r_power, 1);

fprintf("Direct solve L1 difference: %.3e\n", directDifference);
fprintf("Direct solve mass error: %.3e\n\n", directMassError);

comparisonTable = table( ...
    pageNumbers, ...
    r_power, ...
    r_direct, ...
    abs(r_power - r_direct), ...
    'VariableNames', {'Page', 'PowerIteration', 'DirectSolve', 'AbsoluteDifference'} ...
);

disp("Direct solve comparison:");
disp(comparisonTable);

%% Probability Mass Figure

figure;
plot(1:length(massErrors), massErrors, '-o', 'LineWidth', 1.2);
xlabel('Iteration');
ylabel('|sum(r^{(k+1)}) - 1| before renormalization');
title('Probability Mass Error Before Renormalization');
grid on;

saveas(gcf, 'pagerank_mass_error.png');

%% Local Functions

function [S, dangling, outDegree] = build_sparse_transition_matrix(A)

    n = size(A, 1);

    outDegree = full(sum(A, 1));
    dangling = (outDegree == 0);

    invOutDegree = zeros(1, n);
    invOutDegree(~dangling) = 1 ./ outDegree(~dangling);

    Dinv = spdiags(invOutDegree(:), 0, n, n);
    S = A * Dinv;

end

function [r, errors, massErrors, iter] = pagerank_power_sparse(S, dangling, alpha, v, TOL, maxIter)

    r = v;

    errors = zeros(maxIter, 1);
    massErrors = zeros(maxIter, 1);

    for k = 1:maxIter
        danglingMass = sum(r(dangling));

        r_new = alpha * (S * r) ...
              + alpha * danglingMass * v ...
              + (1 - alpha) * v;

        massErrors(k) = abs(sum(r_new) - 1);

        r_new = r_new / sum(r_new);

        err = norm(r_new - r, 1);
        errors(k) = err;

        r = r_new;

        if err < TOL
            iter = k;
            errors = errors(1:k);
            massErrors = massErrors(1:k);
            return;
        end
    end

    iter = maxIter;
    errors = errors(1:maxIter);
    massErrors = massErrors(1:maxIter);

end