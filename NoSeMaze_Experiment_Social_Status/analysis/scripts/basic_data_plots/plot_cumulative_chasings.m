% ========================================================================
% 📌 Script: plot_cumulative_chasings.m
% 📅 Date: 18.03.2025
% 🖊️ Author: Jonathan Reinwald
% 📌 Description:
% - Visualizes cumulative active chasings and times being chased per cohort
% - Computes summary statistics (mean ± SD)
% - Performs statistical comparisons at each rank
% ========================================================================

close all; clearvars;

% ========================================================================
% 🔧 USER SETTINGS
% ========================================================================
% Choose sorting method:
% Option 1: Sort all data based on active chasings (winner-based)
% sorting = 'winner_based';
% Option 2: Sort wins and losses separately
sorting = 'separately';

% Select day range
range_name = 'D1_21';

% ========================================================================
% 📂 LOAD DATA / DEFINE OUTPUT DIRECTORY
% ========================================================================
% Load NoSeMaze dataset (ensure 'summary_data' is available)
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory
summary_file_dir = fullfile(main_dir, 'data', 'processed', 'cross_cohort_files');
load(fullfile(summary_file_dir, 'summary_data.mat'));

% Load cohort information
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));
% Select cohorts that are included in both tube and lickport analyses
cohortsTbl = cohortsTbl((cohortsTbl.use_tube == 1), :);

%% Filter input data (cohort selection)
fields = fieldnames(summary_data.tube);
for field_idx = 1:numel(fields)
    summary_data.tube.(fields{field_idx}) = summary_data.tube.(fields{field_idx})(ismember(summary_data.tube.(fields{field_idx}).cohort,cohortsTbl.cohort),:);
end

% Directory for plotting
plot_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort','basic_data_plots','tube');
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

% Extract relevant variables
cohorts = summary_data.tube.(range_name).cohort;
mouse_IDs = summary_data.tube.(range_name).Mouse_RFID;
fraction_wins = summary_data.tube.(range_name).Fraction_Wins_Chasing;
fraction_losses = summary_data.tube.(range_name).Fraction_Losses_Chasing;
Rank_Chasing = summary_data.tube.(range_name).Rank_Chasing;

% Identify unique cohorts
unique_cohorts = unique(cohorts);

% ========================================================================
% 🎨 FIGURE 1: INDIVIDUAL COHORT CURVES (CUMULATIVE CHASINGS)
% ========================================================================
fig = figure('Visible','on');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.6]);
hold on;
title('Cumulative Chasings Across Cohorts');
xlabel({'Ranked Individuals',['(',sorting,')']},'Interpreter','none');
ylabel('Cumulative Number of Chasings','Interpreter','none');
legend_entries = {};

% Generate red-to-blue colormap
cohort_colors = autumn(length(unique_cohorts));

% Initialize storage for group-level analysis
cumulative_wins_all = [];
cumulative_losses_all = [];
RankChasing_sorted_by_FractionWinsChasing = nan(height(cohortsTbl),10);
RankChasing_sorted_by_FractionLossesChasing = nan(height(cohortsTbl),10);

% ========================================================================
% 📊 LOOP THROUGH COHORTS & PLOT CUMULATIVE SUMS
% ========================================================================
for i = 1:length(unique_cohorts)
    cohort_name = unique_cohorts{i};

    % Filter cohort-specific data
    cohort_idx = strcmp(cohorts, cohort_name);
    cohort_wins = fraction_wins(cohort_idx);
    cohort_losses = fraction_losses(cohort_idx);
    cohort_Rank_Chasing = Rank_Chasing(cohort_idx);

    % Remove NaN values
    valid_idx = ~isnan(cohort_wins) & ~isnan(cohort_losses);
    cohort_wins = cohort_wins(valid_idx);
    cohort_losses = cohort_losses(valid_idx);
    cohort_Rank_Chasing = Rank_Chasing(valid_idx);

    % Skip if insufficient data
    if length(cohort_wins) < 2 || length(cohort_losses) < 2
        fprintf('Skipping cohort %s (not enough data)\n', cohort_name);
        continue;
    end

    % Sorting strategy
    [~, sort_idx_wins] = sort(cohort_wins, 'descend');
    cohort_wins_sorted = cohort_wins(sort_idx_wins);
    if strcmp(sorting, 'separately')
        [~, sort_idx_losses] = sort(cohort_losses, 'descend');
        cohort_losses_sorted = cohort_losses(sort_idx_losses);
    else
        [~, sort_idx_losses] = sort(cohort_losses, 'descend');
        cohort_losses_sorted = cohort_losses(sort_idx_wins);
    end
    % Compute cumulative sums
    cumulative_wins = cumsum(cohort_wins_sorted);
    cumulative_losses = cumsum(cohort_losses_sorted);
    cohort_Rank_Chasing_sortedByFractionWins = cohort_Rank_Chasing(sort_idx_wins);
    cohort_Rank_Chasing_sortedByFractionLosses = cohort_Rank_Chasing(sort_idx_losses);

    % Store data for group-level analysis
    cumulative_wins_all = [cumulative_wins_all; padarray(cumulative_wins(:)', [0, 10 - length(cumulative_wins)], NaN, 'post')];
    cumulative_losses_all = [cumulative_losses_all; padarray(cumulative_losses(:)', [0, 10 - length(cumulative_losses)], NaN, 'post')];
    RankChasing_sorted_by_FractionWinsChasing(i,1:length(cohort_Rank_Chasing_sortedByFractionWins)) = cohort_Rank_Chasing_sortedByFractionWins';
    RankChasing_sorted_by_FractionLossesChasing(i,1:length(cohort_Rank_Chasing_sortedByFractionLosses)) = cohort_Rank_Chasing_sortedByFractionLosses';

    % Plot active chasings (solid) & being chased (dashed)
    plot(1:length(cohort_wins_sorted), cumulative_wins, 'LineWidth', 1, 'Color', cohort_colors(i, :));
    plot(1:length(cohort_losses_sorted), cumulative_losses, '--', 'LineWidth', 1, 'Color', cohort_colors(i, :));

    legend_entries{end+1} = [cohort_name ' (Active)'];
    legend_entries{end+1} = [cohort_name ' (Chased)'];
end

legend(legend_entries, 'Location', 'best','FontSize',5);
hold off;

% save the figure
[~, ~] = docDataSrc(fig,fullfile(plot_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(plot_dir, sprintf('CumulativeChasings_IndividualCohortCurves_%s_%s.pdf', sorting, range_name)));
close(fig);

% ========================================================================
% 🎨 FIGURE 2: SHADED ERROR BARS (GROUP-LEVEL MEAN ± SD)
% ========================================================================
fig = figure('Visible','on');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.6]);
hold on;
title('Mean ± SD of Cumulative Chasings');
xlabel({'Ranked Individuals',['(',sorting,')']},'Interpreter','none');
ylabel('Cumulative Number of Chasings','Interpreter','none');

% Compute mean ± SD
mean_cumulative_wins = nanmean(cumulative_wins_all, 1);
std_cumulative_wins = nanstd(cumulative_wins_all, [], 1);
mean_cumulative_losses = nanmean(cumulative_losses_all, 1);
std_cumulative_losses = nanstd(cumulative_losses_all, [], 1);

% Plot shaded error bars
shadedErrorBar(1:10, mean_cumulative_wins(1:10), std_cumulative_wins(1:10), 'lineProps', {'b-', 'LineWidth', 2});
shadedErrorBar(1:10, mean_cumulative_losses(1:10), std_cumulative_losses(1:10), 'lineProps', {'r--', 'LineWidth', 2});

% Reference line (Equal Distribution)
plot([1,10], [0.1,1], 'k--', 'LineWidth', 1.5);
ylim([0,1]); xlim([0,10]);

legend({'Active Chasings', 'Being Chased'}, 'Location', 'northwest');
hold off;

% ========================================================================
% 📊 STATISTICAL TESTING: ADDING SIGNIFICANCE MARKERS TO FIGURE 2
% ========================================================================
p_values = nan(1, 10); p_values_ttest = nan(1, 10); p_values_perm  = nan(1, 10);
for rank = 1:10
    valid_idx = ~isnan(cumulative_wins_all(:, rank)) & ~isnan(cumulative_losses_all(:, rank));

    if sum(valid_idx) >= 2  % Ensure valid sample size
        if strcmp(sorting, 'separately')
            [~, p_values_ttest(rank)] = ttest2(cumulative_wins_all(valid_idx, rank), cumulative_losses_all(valid_idx, rank));
            p_values(rank) = ranksum(cumulative_wins_all(valid_idx, rank), cumulative_losses_all(valid_idx, rank));
            p_values_perm(rank) = permutationTest(cumulative_wins_all(valid_idx, rank),cumulative_losses_all(valid_idx, rank),10000);
            [observed_diff_2(rank), perm_diffs_2(rank,:), p_values_perm_2(rank)] = permutation_test_unpaired(cumulative_wins_all(valid_idx, rank), cumulative_losses_all(valid_idx, rank), 10000, 'mean');
        elseif strcmp(sorting, 'winner_based')
            [~, p_values_ttest(rank)] = ttest(cumulative_wins_all(valid_idx, rank), cumulative_losses_all(valid_idx, rank));
            p_values(rank) = signrank(cumulative_wins_all(valid_idx, rank), cumulative_losses_all(valid_idx, rank));
            [observed_diff_2(rank), perm_diffs_2(rank,:), p_values_perm_2(rank)] = permutation_test_paired(cumulative_wins_all(valid_idx, rank), cumulative_losses_all(valid_idx, rank), 10000, 'mean');
        end
    end
end


% Define significance markers
sig_levels = [0.05, 0.01, 0.001];
sig_markers = {'*', '**', '***'};

% Plot significance markers above shaded error bars
for rank = 1:10
    if p_values_perm_2(rank) < sig_levels(3)  % p < 0.001
        text(rank, mean_cumulative_wins(rank) + 0.05, sig_markers{3}, 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    elseif p_values_perm_2(rank) < sig_levels(2)  % p < 0.01
        text(rank, mean_cumulative_wins(rank) + 0.05, sig_markers{2}, 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    elseif p_values_perm_2(rank) < sig_levels(1)  % p < 0.05
        text(rank, mean_cumulative_wins(rank) + 0.05, sig_markers{1}, 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
end

% ========================================================================
% ✅ DISPLAY STATISTICAL RESULTS
% ========================================================================

% Display results
disp('Paired t-test p-values:');
disp(array2table(p_values_ttest, 'VariableNames', strcat("Rank_", string(1:10))));
disp('Wilcoxon signed-rank or Mann-Whitney U test p-values:');
disp(array2table(p_values, 'VariableNames', strcat("Rank_", string(1:10))));
if strcmp(sorting, 'separately')
    disp('Permutation test p-values:');
    disp(array2table(p_values_perm_2, 'VariableNames', strcat("Rank_", string(1:10))));
end

% Source Data
SourceDataTable = array2table(cumulative_wins_all);
writetable(SourceDataTable,fullfile(plot_dir, ['SourceData_CumulativeWins_AllCohorts.csv']));

% Source Data
SourceDataTable = array2table(cumulative_losses_all);
writetable(SourceDataTable,fullfile(plot_dir, ['SourceData_CumulativeLosses_AllCohorts.csv']));

% save the figure
[~, ~] = docDataSrc(fig,fullfile(plot_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(plot_dir, sprintf('CumulativeChasings_ShadedErrorBars_%s_%s.pdf', sorting, range_name)));
close(fig);

% ========================================================================
% ✅ END OF SCRIPT
% ========================================================================