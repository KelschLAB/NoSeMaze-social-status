% ========================================================================
% Script: plot_compare_distributions_social_traits.m
% Date: 18.03.2025
% Author: Jonathan Reinwald
% Description:
%   This script analyzes the differences between "win" and "loss"
%   fractions in selected behaviors (e.g., Competition, Chasing)
%   using histograms with KDE fits, CDF plots, and the KS test.
%
%   It also saves figures and test results for reproducibility.
% ========================================================================

close all; clearvars;

% ========================================================================
% USER SETTINGS
% ========================================================================
selection = {'Competition', 'Chasing'};
range_name = {'D1_21','D1_End'};

% ========================================================================
% LOAD DATA / DEFINE OUTPUT DIRECTORY
% ========================================================================
% Load NoSeMaze dataset (ensure 'summary_data' is available)
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory
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

% Initialize results table
ks_results = table('Size', [length(selection), 3], ...
    'VariableTypes', {'string', 'double', 'logical'}, ...
    'VariableNames', {'Behavior', 'pValue', 'IsDifferent'});

% ========================================================================
% Starting loop
% ========================================================================
% Over ranges
for rx = 1:length(range_name)

    current_range = range_name{rx};
    % Over selections
    for ix = 1:length(selection)

        behavior = selection{ix};

        % Extract and clean data
        wins = summary_data.tube.(current_range).(['Fraction_Wins_' behavior]);
        losses = summary_data.tube.(current_range).(['Fraction_Losses_' behavior]);
        wins = wins(~isnan(wins));
        losses = losses(~isnan(losses));

        % --- Save source data as CSV ---
        if ~exist(plot_dir, 'dir')
            mkdir(plot_dir);
        end
        % Create a table and save
        T = table(wins(:), losses(:), ...
            'VariableNames', {'Fraction_Wins', 'Fraction_Losses'});
        filename = fullfile(plot_dir, sprintf('SourceData_%s_%s.csv', behavior, current_range));
        writetable(T, filename);
        disp(['Saved source data to: ' filename]);

        % KS test
        [h, p] = kstest2(wins, losses);
        ks_results.Behavior(ix) = behavior;
        ks_results.pValue(ix) = p;
        ks_results.IsDifferent(ix) = h;

        % Prepare figure
        fig = figure('Name', ['Distribution Comparison: ' behavior], 'Color', 'w');

        % Shared range
        all_data = [wins; losses];
        x_min = min(all_data);
        x_max = max(all_data);
        x_range = linspace(x_min, x_max, 200);

        % --- Subplot 1: Histogram + KDE ---
        subplot(1,2,1);
        histogram(wins, 'BinWidth', 0.02, ...%'Normalization', 'pdf', ...
            'EdgeColor', 'none', 'FaceColor', [0.2 0.4 1], 'FaceAlpha', 0.4);
        hold on;
        histogram(losses, 'BinWidth', 0.02, ...%'Normalization', 'pdf', ...
            'EdgeColor', 'none', 'FaceColor', [1 0.2 0.2], 'FaceAlpha', 0.4);
        % Bin width (same as histogram)
        binWidth = 0.02;

        % Scale KDE to match total count
        [f_wins, x_wins] = ksdensity(wins);
        f_wins_scaled = f_wins * length(wins) * binWidth;
        plot(x_wins, f_wins_scaled, 'b-', 'LineWidth', 2);

        [f_losses, x_losses] = ksdensity(losses);
        f_losses_scaled = f_losses * length(losses) * binWidth;
        plot(x_losses, f_losses_scaled, 'r-', 'LineWidth', 2);
        xlabel('Fraction');
        ylabel('Count');
        title(['PDF + KDE: ' behavior]);
        legend({'Wins Histogram', 'Losses Histogram', 'Wins Fit (scaled)', 'Losses Fit (scaled)'});
        xlim([0 0.6]);
        box off;
        axis square;


        % Format test result text
        sig_text = 'No';
        if h == 1
            sig_text = 'Yes';
        end
        ks_text = sprintf('KS Test:\np = %.4f\nSignificant: %s', p, sig_text);

        % Add annotation to top right of PDF subplot
        yl = ylim;
        text(x_max - 0.05 * range([x_min x_max]), yl(2)*0.9, ks_text, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'white', ...
            'EdgeColor', [0.8 0.8 0.8]);

        grid on;

        % --- Subplot 2: Empirical CDFs ---
        subplot(1,2,2);
        cdfplot(wins); hold on;
        cdfplot(losses);
        xlabel('Fraction');
        ylabel('Cumulative Probability');
        title(['Empirical CDFs: ' behavior]);
        legend('Wins', 'Losses');
        xlim([0 0.6]);
        grid on;
        box off;
        axis square;

        % --- Save figure ---
        [~, ~] = docDataSrc(fig,fullfile(plot_dir),mfilename('fullpath'),logical(1));
        exportgraphics(fig, fullfile(plot_dir, sprintf('DistComp_%s_%s.pdf', behavior,current_range)));
        disp(['Saved figure to: ' fullfile(plot_dir, sprintf('DistComp_%s_%s.pdf', behavior,current_range))]);
        close(fig);
    end

    % Save KS results
    results_file = fullfile(plot_dir, sprintf('KS_test_results_%s_.csv',current_range));
    writetable(ks_results, results_file);
    disp(['Saved KS test results to: ' results_file]);
end
