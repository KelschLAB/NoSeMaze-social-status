function [observed_diff, perm_diffs, p_value] = permutation_test_unpaired(vec1, vec2, n_perm, method)
    % Perform a permutation test for two unpaired vectors (Mean or Median-based)
    %
    % Inputs:
    %   vec1   - First sample (vector)
    %   vec2   - Second sample (vector)
    %   n_perm - Number of permutations
    %   method - 'mean' (default) or 'median' for comparison
    %
    % Outputs:
    %   observed_diff - Observed difference in means or medians
    %   perm_diffs    - Vector of all permuted differences
    %   p_value       - Computed p-value

    if nargin < 4
        method = 'mean';  % Default method is 'mean'
    end

    % Compute the observed difference
    switch method
        case 'mean'
            observed_diff = mean(vec1) - mean(vec2);
        case 'median'
            observed_diff = median(vec1) - median(vec2);
        otherwise
            error('Invalid method. Choose ''mean'' or ''median''.');
    end

    % Pool both vectors together
    combined = [vec1; vec2];
    len_vec1 = length(vec1);
    
    % Initialize array to store permuted differences
    perm_diffs = zeros(n_perm, 1);
    
    for i = 1:n_perm
        % Randomly permute the combined data
        shuffled = combined(randperm(length(combined)));

        % Split into new "random" groups
        new_vec1 = shuffled(1:len_vec1);
        new_vec2 = shuffled(len_vec1+1:end);

        % Compute the new difference based on the selected method
        switch method
            case 'mean'
                perm_diffs(i) = mean(new_vec1) - mean(new_vec2);
            case 'median'
                perm_diffs(i) = median(new_vec1) - median(new_vec2);
        end
    end

    % Compute p-value: Two-tailed test (how often is permuted difference as extreme as observed?)
    p_value = mean(abs(perm_diffs) >= abs(observed_diff));

    % Display results
    fprintf('Permutation Test (%s-based):\n', method);
    fprintf('Observed Difference: %.4f\n', observed_diff);
    fprintf('Permutation Test p-value: %.4f (with %d permutations)\n', p_value, n_perm);
end
