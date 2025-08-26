function [observed_diff, perm_diffs, p_value] = permutation_test_paired(vec1, vec2, n_perm, method)
% Perform a permutation test for paired data (mean or median-based)
%
% Inputs:
%   vec1   - First condition (e.g., pre-test)
%   vec2   - Second condition (e.g., post-test)
%   n_perm - Number of permutations
%   method - 'mean' (default) or 'median'
%
% Outputs:
%   observed_diff - Observed difference (mean or median of paired diffs)
%   perm_diffs    - All permuted differences
%   p_value       - Two-tailed p-value

    if nargin < 4
        method = 'mean';
    end

    if length(vec1) ~= length(vec2)
        error('Input vectors must have the same length for paired test.');
    end

    % Compute the paired differences
    diffs = vec1 - vec2;

    % Observed test statistic
    switch method
        case 'mean'
            observed_diff = mean(diffs);
        case 'median'
            observed_diff = median(diffs);
        otherwise
            error('Invalid method. Choose ''mean'' or ''median''.');
    end

    % Initialize permutation differences
    perm_diffs = zeros(n_perm, 1);

    for i = 1:n_perm
        % Randomly flip signs
        signs = randi([0 1], size(diffs)) * 2 - 1;  % Generates +1 or -1
        flipped = diffs .* signs;

        % Calculate permuted test statistic
        switch method
            case 'mean'
                perm_diffs(i) = mean(flipped);
            case 'median'
                perm_diffs(i) = median(flipped);
        end
    end

    % Two-tailed p-value
    p_value = mean(abs(perm_diffs) >= abs(observed_diff));

    % Display results
    fprintf('Paired Permutation Test (%s-based):\n', method);
    fprintf('Observed Difference: %.4f\n', observed_diff);
    fprintf('Permutation Test p-value: %.4f (with %d permutations)\n', p_value, n_perm);
end
