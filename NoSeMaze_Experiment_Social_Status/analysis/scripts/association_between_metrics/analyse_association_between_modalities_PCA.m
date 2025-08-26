%% analyze_association_between_modalities_PCA.m
% Author: Jonathan Reinwald
% Date: April 2025

% Purpose:
% Analyze associations between social and cognitive features with a PCA.

%% Pre-Clearing
clear; clc; close all;

%% Set Directories
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory
processed_dir = fullfile(main_dir, 'data', 'processed');
results_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort', 'associations_between_modalities');
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');

% Add required paths
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));
addpath(genpath(fullfile(main_dir, 'src', 'analysis', 'associations_between_modalities')));
addpath(genpath(fullfile(main_dir, 'src', 'analysis', 'associations_between_metrics')));

%% Load Summary Data
if ~isfile(summary_file)
    error('Summary file not found. Please run `main_build_cross_cohort_summary` first.');
end
load(summary_file, 'summary_data');

%% Cohort Selection
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));

%% Define Parameters
data_types = {'tube','lickport'}; % Analyze tube metrics
X_metrics.tube = {'DSz_Competition', 'Fraction_Wins_Chasing', 'Fraction_Losses_Chasing'};
Y_metrics.tube = {'DSz_Competition', 'Fraction_Wins_Chasing', 'Fraction_Losses_Chasing'};
X_metrics.lickport = {'correct_hit_rate_allTrials','correct_rejection_rate_allTrials',...
    'baseline_rate_mean_omitfirst','cs_plus_modulation_peak',...
    'cs_plus_switch_latency_at_cs_rev1','cs_minus_switch_latency_at_cs_rev1','cs_plus_switch_latency_at_cs_rev2','cs_minus_switch_latency_at_cs_rev2',...
    'cs_plus_switch_latency_at_cs_rev3','cs_minus_switch_latency_at_cs_rev3','cs_plus_switch_latency_at_cs_rev4','cs_minus_switch_latency_at_cs_rev4','cs_plus_switch_latency_at_cs_median','cs_minus_switch_latency_at_cs_median'};
Y_metrics.lickport = {'correct_hit_rate_allTrials','correct_rejection_rate_allTrials',...
    'baseline_rate_mean_omitfirst','cs_plus_modulation_peak',...
    'cs_plus_switch_latency_at_cs_rev1','cs_minus_switch_latency_at_cs_rev1','cs_plus_switch_latency_at_cs_rev2','cs_minus_switch_latency_at_cs_rev2',...
    'cs_plus_switch_latency_at_cs_rev3','cs_minus_switch_latency_at_cs_rev3','cs_plus_switch_latency_at_cs_rev4','cs_minus_switch_latency_at_cs_rev4','cs_plus_switch_latency_at_cs_median','cs_minus_switch_latency_at_cs_median'};
day_ranges = {'D1_End'}; % Select day_ranges of interest
LME_covariates = {'genotype'}; % Optional covariates: If none: LME_covariates = {''}


% Output directory
plot_output_dir = fullfile(results_dir, 'overview_PCA');
if ~isfolder(plot_output_dir)
    mkdir(plot_output_dir);
end

%% Iterate Through Ray Ranges, Data Types, ...
% Loop over day ranges
for range_idx = 1:length(day_ranges)
    range_name = day_ranges{range_idx};

    % Combined metrics
    X_metrics_combined.(range_name) = [];
    Y_metrics_combined.(range_name) = [];
    data_table_combined.(range_name) = [];

    % Loop over data type (tube and lickport data)
    for type_idx = 1:length(data_types)
        data_type = data_types{type_idx};

        % select cohorts
        cohortsTbl = cohortsTbl((cohortsTbl.use_tube == 1) & (cohortsTbl.use_lickport == 1), :);

        % Filter input data (cohort selection)
        fields = fieldnames(summary_data.(data_type));
        for field_idx = 1:numel(fields)
            summary_data.(data_type).(fields{field_idx}) = summary_data.(data_type).(fields{field_idx})(ismember(summary_data.(data_type).(fields{field_idx}).cohort,cohortsTbl.cohort),:);
        end

        % Combined metrics
        X_metrics_combined.(range_name) = [X_metrics_combined.(range_name),X_metrics.(data_type)];
        Y_metrics_combined.(range_name) = [Y_metrics_combined.(range_name),Y_metrics.(data_type)];
        if isempty(data_table_combined.(range_name))
            data_table_combined.(range_name) = summary_data.(data_type).(range_name)(:,['Mouse_RFID','cohort','genotype',X_metrics.(data_type)]);
        else
            data_table_combined.(range_name) = join(data_table_combined.(range_name),summary_data.(data_type).(range_name)(:,['Mouse_RFID','cohort',X_metrics.(data_type)]),'Keys',{'Mouse_RFID','cohort'});
        end
    end

    % Define datTable
    dataTable = data_table_combined.(range_name);

    % Extract numeric data matrix and variable names
    dataToNormalize = dataTable{:, X_metrics_combined.(range_name)};
    varNames = dataTable.Properties.VariableNames(X_metrics_combined.(range_name));

    % Preallocate transformed data
    transformedData = zeros(size(dataToNormalize));

    % Plot histograms toggle
    plotHistograms = true;

    % Step 1: Adaptive transformation based on skewness
    epsilon = 0.01;
    nVars = size(dataToNormalize, 2);

    % Histogram layout
    plotsPerFigure = 10;  % Total subplots per figure
    subplotsUsed = 0;    % Counter

    if plotHistograms
        figCount = 1;
        fig=figure('Name', ['Histograms Set ' num2str(figCount)], 'Position', [20 20 800 1000]);
    end

    for i = 1:nVars
        x = dataToNormalize(:, i);
        x_nonan = x(~isnan(x));
        s_pre = skewness(x_nonan);

        % Default: no transform
        x_trans = x;
        s_post = s_pre;

        if abs(s_pre) > 1
            shift = abs(min(x_nonan)) + epsilon;
            if s_pre > 1
                x_log = log(x + shift);
                x_sqrt = sqrt(x + shift);
            else
                reflected = max(x_nonan) + shift - x;
                x_log = log(reflected);
                x_sqrt = sqrt(reflected);
            end

            s_log = skewness(x_log(~isnan(x_log)));
            s_sqrt = skewness(x_sqrt(~isnan(x_sqrt)));

            if abs(s_log) < abs(s_sqrt)
                x_trans = x_log;
                s_post = s_log;
            else
                x_trans = x_sqrt;
                s_post = s_sqrt;
            end
        end

        transformedData(:, i) = x_trans;

        if plotHistograms
            % Plot before
            subplot(5, 2, subplotsUsed + 1);
            histogram(x_nonan, 20, 'EdgeColor', 'none');
            title({['Before: ' varNames{i}], ['Skew = ' num2str(s_pre, '%.2f')]}, 'Interpreter', 'none');
            xlabel(varNames{i}, 'Interpreter', 'none');
            ylabel('Count', 'Interpreter', 'none');
            subplotsUsed = subplotsUsed + 1;

            % Check if new figure is needed
            if subplotsUsed >= plotsPerFigure
                % Save current figure
                figFilename = sprintf('Histograms_%s_Set%d.pdf', range_name, figCount);
                figPath = fullfile(plot_output_dir, figFilename);
                [annot, srcInfo] = docDataSrc(fig, plot_output_dir, mfilename('fullpath'), true);
                exportgraphics(fig, figPath);

                % Open new figure
                figCount = figCount + 1;
                fig = figure('Name', ['Histograms Set ' num2str(figCount)], 'Position', [20 20 800 1000]);
                subplotsUsed = 0;
            end

            % Plot after
            subplot(5, 2, subplotsUsed + 1);
            histogram(x_trans(~isnan(x_trans)), 20, 'EdgeColor', 'none');
            title({['After: ' varNames{i}], ['Skew = ' num2str(s_post, '%.2f')]}, 'Interpreter', 'none');
            xlabel(varNames{i}, 'Interpreter', 'none');
            ylabel('Count', 'Interpreter', 'none');
            subplotsUsed = subplotsUsed + 1;

            % Save and open new figure if full
            if subplotsUsed >= plotsPerFigure && i < nVars
                figFilename = sprintf('Histograms_%s_Set%d.pdf', range_name, figCount);
                figPath = fullfile(plot_output_dir, figFilename);
                [annot, srcInfo] = docDataSrc(fig, plot_output_dir, mfilename('fullpath'), true);
                exportgraphics(fig, figPath);

                figCount = figCount + 1;
                fig = figure('Name', ['Histograms Set ' num2str(figCount)], 'Position', [20 20 800 1000]);
                subplotsUsed = 0;
            end
        end
    end

    % Save last figure if not already saved
    if plotHistograms && subplotsUsed > 0
        figFilename = sprintf('Histograms_%s_Set%d.pdf', range_name, figCount);
        figPath = fullfile(plot_output_dir, figFilename);
        [annot, srcInfo] = docDataSrc(fig, plot_output_dir, mfilename('fullpath'), true);
        exportgraphics(fig, figPath);
    end

    % Step 2: Z-score standardization
    zscoredData = (transformedData - mean(transformedData, 1, 'omitnan')) ./ std(transformedData, 0, 1, 'omitnan');

    % Step 3: Imputation
    for i = 1:size(zscoredData, 2)
        x = zscoredData(:, i);
        x(isnan(x)) = mean(x, 'omitnan');
        zscoredData(:, i) = x;
    end

    % Perform PCA
    [coeff, score, latent, tsquared, explained, mu] = pca(zscoredData);

    % === Plot 1: Explained Variance ===
    fig = figure;
    pareto(explained);
    xlabel('Principal Component', 'Interpreter', 'none');
    ylabel('Variance Explained (%)', 'Interpreter', 'none');
    title('Variance Explained by Principal Components', 'Interpreter', 'none');

    filename = sprintf('PCA_VarianceExplained_%s.pdf', range_name);
    [annot, srcInfo] = docDataSrc(fig, plot_output_dir, mfilename('fullpath'), true);
    exportgraphics(fig, fullfile(plot_output_dir, filename));

    % --- Define domain groups
    socialVars = {'DSz_Competition', 'Fraction_Wins_Chasing', 'Fraction_Losses_Chasing'};
    cognitiveVars = setdiff(varNames, socialVars);

    % Logical indices for column selection
    socialIdx = ismember(varNames, socialVars);
    cognitiveIdx = ismember(varNames, cognitiveVars);

    % --- Permutation test for domain loading difference significance (PCs 1–5) ---
    nPerm = 10000;
    p_perm_vec = zeros(1, 5);  % Store p-values for later annotation

    for pc = 1:5
        observedDiff = mean(abs(coeff(socialIdx, pc))) - mean(abs(coeff(cognitiveIdx, pc)));
        combined = [abs(coeff(socialIdx, pc)); abs(coeff(cognitiveIdx, pc))];
        n1 = length(socialIdx(socialIdx==1));
        permDiffs = zeros(nPerm,1);
        for i = 1:nPerm
            idx = randperm(length(combined));
            permDiffs(i) = mean(combined(idx(1:n1))) - mean(combined(idx(n1+1:end)));
        end
        p_perm = mean(abs(permDiffs) >= abs(observedDiff));
        p_perm_vec(pc) = p_perm;  % Store
        fprintf('Permutation test PC%d: Observed diff = %.3f, p = %.4f\n', pc, observedDiff, p_perm);
    end

    % === Plot 2: PCA Loadings ===
    fig = figure;
    bar(coeff(:, 1:3), 'grouped', 'EdgeColor', 'none');
    legend({'PC1','PC2','PC3'}, 'Location', 'northeast', 'Interpreter', 'none');
    xticks(1:length(varNames));
    xticklabels(varNames);
    xtickangle(90);
    ylabel('Loading Value', 'Interpreter', 'none');
    title('Loadings of First 3 Principal Components', 'Interpreter', 'none');

    % Disable LaTeX interpretation for tick labels
    ax = gca;
    ax.TickLabelInterpreter = 'none';
    ax.XTickLabelRotation = 45;

    % Add p-values to plot (as annotations above the legend)
    text(0.05, 1.1, sprintf('p_{PC1} = %.4f', p_perm_vec(1)), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'Interpreter', 'none');
    text(0.05, 1.05, sprintf('p_{PC2} = %.4f', p_perm_vec(2)), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'Interpreter', 'none');
    text(0.05, 1.00, sprintf('p_{PC3} = %.4f', p_perm_vec(3)), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'Interpreter', 'none');

    % Source data
    writetable(table(varNames',coeff), fullfile(plot_output_dir, 'PCA_coeff.csv'));

    % Save figure
    filename = sprintf('PCA_Loadings_%s.pdf', range_name);
    [annot, srcInfo] = docDataSrc(fig, plot_output_dir, mfilename('fullpath'), true);
    exportgraphics(fig, fullfile(plot_output_dir, filename));


    % === Show top contributors (console output only) ===
    for pc = 1:5
        [~, idx] = sort(abs(coeff(:, pc)), 'descend');
        disp(['Top variables contributing to PC' num2str(pc) ':']);
        disp(varNames(idx(1:5)));
    end

    % % === Plot 3: Correlation Matrix ===
    % fig = figure;
    % imagesc(corr(zscoredData(:, socialIdx), zscoredData(:, cognitiveIdx),'type','Spearman')); % CAVE: also NAN rows excluded where data is complete!
    % xticks(1:length(cognitiveVars));
    % yticks(1:length(socialVars));
    % xticklabels(cognitiveVars);
    % yticklabels(socialVars);
    % title('Cross-domain Correlations (Social vs. Cognitive)', 'Interpreter', 'none');
    % xlabel('Cognitive Traits', 'Interpreter', 'none');
    % ylabel('Social Traits', 'Interpreter', 'none');
    % % colorbar
    % ax=gca; 
    % crameri('vik'); %ax.Colormap=redbluecmap;
    % ax.CLim=[-1,1];
    % colorbar;
    % ax.TickLabelInterpreter = 'none';
    % 
    % filename = sprintf('CorrelationMatrix_%s.pdf', range_name);
    % [annot, srcInfo] = docDataSrc(fig, plot_output_dir, mfilename('fullpath'), true);
    % exportgraphics(fig, fullfile(plot_output_dir, filename));

    % Multivariate regression
    X = zscoredData(:, socialIdx);
    Y = zscoredData(:, cognitiveIdx);
    B = X \ Y;
    Y_pred = X * B;

    % R² output
    R2 = 1 - sum((Y - Y_pred).^2) ./ sum((Y - mean(Y)).^2);
    for i = 1:length(cognitiveVars)
        fprintf('R^2 for %s: %.4f\n', cognitiveVars{i}, R2(i));
    end

    % === Plot 4: PCA Clusters ===
    numClusters = 3;
    clusterLabels = kmeans(zscoredData, numClusters);
    dataTable.Cluster = clusterLabels;

    fig = figure;
    gscatter(score(:,1), score(:,2), clusterLabels);
    xlabel('PC1', 'Interpreter', 'none');
    ylabel('PC2', 'Interpreter', 'none');
    title('Behavioral Clusters in PCA Space', 'Interpreter', 'none');
    legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Interpreter', 'none');

    % Disable LaTeX interpretation for tick labels
    ax = gca;
    ax.TickLabelInterpreter = 'none';

    filename = sprintf('PCA_Clusters_%s.pdf', range_name);
    [annot, srcInfo] = docDataSrc(fig, plot_output_dir, mfilename('fullpath'), true);
    exportgraphics(fig, fullfile(plot_output_dir, filename));

    % === Save Cluster Profiles (Table + Heatmap) ===
    % High dominance, low flexibility
    DSz_idx = find(strcmp(varNames, 'DSz_Competition'));
    PC1_scores = score(:, 1);
    highDomLowCog = dataTable(zscoredData(:, DSz_idx) > 1 & PC1_scores < -1, :);
    disp('Mismatched profiles: high dominance, low cognitive flexibility');
    disp(highDomLowCog);

    % Low dominance, high flexibility
    lowDomHighCog = dataTable(zscoredData(:, DSz_idx) < -1 & PC1_scores > 1, :);
    disp('Mismatched profiles: low dominance, high cognitive flexibility');
    disp(lowDomHighCog);

    % Cluster summary
    clusterSummary = [];
    for c = 1:numClusters
        clusterMeans = mean(zscoredData(clusterLabels == c, :));
        clusterSummary = [clusterSummary; clusterMeans];
        fprintf('\nCluster %d mean trait values:\n', c);
        disp(array2table(clusterMeans, 'VariableNames', varNames));
    end

    % Save cluster means as table
    clusterSummaryTable = array2table(clusterSummary, 'VariableNames', varNames);
    clusterSummaryTable.Cluster = (1:numClusters)';
    writetable(clusterSummaryTable, fullfile(plot_output_dir, sprintf('ClusterSummary_%s.csv', range_name)));

    % === Plot 5: Cluster Heatmap ===
    fig = figure;
    imagesc(clusterSummary);
    colorbar;
    colormap('parula');
    title('Cluster Mean Profiles', 'Interpreter', 'none');
    xlabel('Traits', 'Interpreter', 'none');
    ylabel('Cluster', 'Interpreter', 'none');
    xticks(1:length(varNames));
    xticklabels(varNames);
    xtickangle(90);
    yticks(1:numClusters);
    yticklabels(arrayfun(@(x) ['Cluster ' num2str(x)], 1:numClusters, 'UniformOutput', false));

    % Disable LaTeX interpretation for tick labels
    ax = gca;
    ax.TickLabelInterpreter = 'none';

    filename = sprintf('ClusterHeatmap_%s.pdf', range_name);
    [annot, srcInfo] = docDataSrc(fig, plot_output_dir, mfilename('fullpath'), true);
    exportgraphics(fig, fullfile(plot_output_dir, filename));



end
