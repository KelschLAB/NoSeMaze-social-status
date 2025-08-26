function [p, observed_stat, perm_stats] = permutation_rank_test(group1, group2, n_perms, paired)
    % permutation_rank_test: Performs a permutation rank test (paired or unpaired)
    % Inputs:
    %   - group1, group2: Numeric vectors of equal length (paired) or any length (unpaired)
    %   - n_perms: Number of permutations
    %   - paired: Boolean flag (true for paired, false for unpaired)
    % Outputs:
    %   - p: Two-sided p-value
    %   - observed_stat: Observed test statistic
    %   - perm_stats: Permutation distribution of test statistics

    if paired
        % Ensure equal length
        assert(length(group1) == length(group2), 'Groups must be the same length for paired test.');
        
        % Compute differences and ranks
        diffs = group1 - group2;
        non_zero_diffs = diffs(diffs ~= 0); % Wilcoxon ignores zeros
        ranks = tiedrank(abs(non_zero_diffs));
        signs = sign(non_zero_diffs);

        % Observed Wilcoxon signed-rank statistic
        observed_stat = sum(ranks .* signs);

        % Permutation: randomly flip signs
        perm_stats = zeros(n_perms, 1);
        for i = 1:n_perms
            flip_signs = sign(rand(length(non_zero_diffs), 1) - 0.5); % random +/-1
            perm_stats(i) = sum(ranks .* flip_signs);
        end
    else
        % Unpaired test (Mann-Whitney U test logic)
        combined = [group1(:); group2(:)];
        n1 = length(group1);
        n2 = length(group2);

        ranks = tiedrank(combined);
        observed_stat = sum(ranks(1:n1));  % Rank sum for group1

        perm_stats = zeros(n_perms, 1);
        for i = 1:n_perms
            perm_idx = randperm(n1 + n2);
            perm_ranks = ranks(perm_idx);
            perm_stats(i) = sum(perm_ranks(1:n1));
        end
    end

    % Two-sided p-value
    p = mean(abs(perm_stats - mean(perm_stats)) >= abs(observed_stat - mean(perm_stats)));
end
