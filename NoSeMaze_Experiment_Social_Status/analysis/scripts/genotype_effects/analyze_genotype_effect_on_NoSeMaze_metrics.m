%% analyze_genotype_effect_on_NoSeMaze_metrics.m
% Author: Jonathan Reinwald
% Date: April 2025

% Purpose:
% Compare genotype (WT vs OxtKO) effects across all NoSeMaze metrics in both tube and lickport datasets.

%% Pre-Clearing
clear; clc; close all;

%% Set Directories
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory
processed_dir = fullfile(main_dir, 'data', 'processed');
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');

% Add helpers if needed
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));

%% Load Summary Data
if ~isfile(summary_file)
    error('Summary file not found. Please run `main_build_cross_cohort_summary` first.');
end
load(summary_file, 'summary_data');

%% Define Parameters
data_types = {'tube', 'lickport'}; % Now doing both
day_ranges.tube = {'D1_7','D8_14','D15_21','D1_21','D1_14'};
day_ranges.lickport = {'D1_End'};

% Define which metrics to analyze per data type
NoSeMaze_metrics.tube = {'DSz_Competition','Fraction_Wins_Chasing', 'Fraction_Losses_Chasing','N_Events_Chasing','detection_count','detection_percentage','detections_per_day'};
%     'Fraction_Wins_Competition', 'Fraction_Losses_Competition','N_Wins_Competition', 'N_Losses_Competition', 'N_Events_Competition','DSz_Chasing','Fraction_Wins_Chasing', 'Fraction_Losses_Chasing', 'N_Wins_Chasing', 'N_Losses_Chasing','N_Events_Chasing',...
%     'detection_count','detection_percentage','detections_per_day'};
NoSeMaze_metrics.lickport = {'correct_hit_rate_allTrials','correct_rejection_rate_allTrials', ...
    'baseline_rate_mean_omitfirst','cs_plus_modulation_peak', ...
    'cs_plus_switch_latency_at_cs_median','cs_minus_switch_latency_at_cs_median'};
    % 'cs_plus_switch_latency_at_cs_rev1','cs_minus_switch_latency_at_cs_rev1', ...
    % 'cs_plus_switch_latency_at_cs_rev2','cs_minus_switch_latency_at_cs_rev2', ...
    % 'cs_plus_switch_latency_at_cs_rev3','cs_minus_switch_latency_at_cs_rev3', ...
    % 'cs_plus_switch_latency_at_cs_rev4','cs_minus_switch_latency_at_cs_rev4', ...


genotype_field = 'genotype'; % field name for genotype

 % Define result directory
results_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort', 'genotype_effects',genotype_field);

%% Load cohort information and select valid cohorts
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));

%% Plot colors
plot_colors = [113 189 134; 178 170 209; 218 147 74] / 255;  % Normalize to [0,1]

%% Method
% methods = {'mean','median','ranks'};
methods = {'median'};


%% Loop through methods
for mx = 1:length(methods)
    myMethod = methods{mx};
    %% Loop through data types
    for dt = 2:length(data_types)
        data_type = data_types{dt};

        % Select valid cohorts for each data type
        if strcmp(data_type, 'tube')
            valid_cohorts = cohortsTbl.cohort(cohortsTbl.use_tube == 1);
        elseif strcmp(data_type, 'lickport')
            valid_cohorts = cohortsTbl.cohort(cohortsTbl.use_lickport == 1);
        else
            error('Unknown data type: %s', data_type);
        end

        % Now loop through day ranges for this data type
        for day_idx = 1:numel(day_ranges.(data_type))
            day_range = day_ranges.(data_type){day_idx};

            % Check if day range exists
            if ~isfield(summary_data.(data_type), day_range)
                warning('Day range %s not found in %s data. Skipping...', day_range, data_type);
                continue;
            end

            % Extract data
            tbl = summary_data.(data_type).(day_range);

            % Filter valid cohorts
            tbl = tbl(ismember(tbl.cohort, valid_cohorts), :);

            % Filter valid genotypes
            valid_idx = ismember(tbl.(genotype_field), {'WT', 'OxtKO'});
            tbl = tbl(valid_idx,:);

            % Prepare output
            output_dir = fullfile(results_dir, data_type, day_range);
            if ~isfolder(output_dir)
                mkdir(output_dir);
            end

            % Initialize results table
            Results = table();

            %% Loop through each metric
            for metric_idx = 1:length(NoSeMaze_metrics.(data_type))
                metric_name = NoSeMaze_metrics.(data_type){metric_idx};

                % Check if metric exists in table
                if ~ismember(metric_name, tbl.Properties.VariableNames)
                    warning('Metric %s not found in %s data. Skipping...', metric_name, data_type);
                    continue;
                end

                % Extract data
                yData = tbl.(metric_name);
                group = categorical(tbl.(genotype_field)); % WT = 1, KO = 2
                animal_ID = tbl.Mouse_RFID;

                % Remove missing data
                valid_entries = ~isnan(yData);
                yData = yData(valid_entries);
                group = group(valid_entries);
                animal_ID = animal_ID(valid_entries);

                %% Prepare Data
                yData = yData(:);
                group = group(:);
                animal_ID = animal_ID(:);

                % Fix animal_ID type
                animal_ID = string(animal_ID);

                group_numeric = double(group == 'OxtKO'); % WT=0, KO=1

                % Run full analysis (Permutation + lme)
                [AnalysisResults, lme_model] = permutation_and_lme_analysis(yData, group, animal_ID, 'n_perm', 10000, 'method', myMethod, 'run_lme', true);

                %% LME model check
                if strcmp(myMethod,'mean')
                    check_LME_assumptions(lme_model, ...
                        'fig_name', sprintf('%s_by_Genotype_%s_%s.pdf', metric_name, day_range, myMethod), ...
                        'plot_save_path', output_dir,...
                        'run_robust', true);
                end

                % Choose correct effect size label
                switch myMethod
                    case 'mean'
                        effect_label = sprintf("Cohen''s d = %.2f", AnalysisResults.cohen_d);
                    case {'median', 'ranks'}
                        effect_label = sprintf("Cliff''s delta = %.2f", AnalysisResults.cliffs_delta);
                    otherwise
                        effect_label = '';
                end

                % Build plot title with effect size
                plot_title = sprintf('%s by Genotype (%s)\nDiffObs p=%.3f | Perm p=%.3f | LME p=%.3f | LMEBoxCox p=%.3f | %s', ...
                    metric_name, day_range, ...
                    AnalysisResults.diff_obs, ...
                    AnalysisResults.pval_permutation, ...
                    AnalysisResults.pval_lme, ...
                    AnalysisResults.pval_lme_log, ...
                    effect_label);

                %% Plot
                fig = figure('visible', 'off');
                subplot(2,2,3)
                plot_data = yData;
                grouping_data = double(group == 'OxtKO') + 1; % 1 = WT, 2 = KO

                % Plot boxplot
                h = boxplot(plot_data, grouping_data, ...
                    'Symbol', 'o', 'OutlierSize', 0.01, 'Whisker', 1.5, ...
                    'colors', plot_colors, 'boxstyle', 'outline');

                % Customize boxplot elements
                set(h(1,:), 'LineWidth', 1.5, 'LineStyle', '-'); % Upper Whisker
                set(h(2,:), 'LineWidth', 1.5, 'LineStyle', '-'); % Lower Whisker
                set(h(3,:), 'LineWidth', 1.5, 'LineStyle', '-'); % Upper Adjacent Value
                set(h(4,:), 'LineWidth', 1.5, 'LineStyle', '-'); % Lower Adjacent Value
                set(h(5,:), 'LineWidth', 1.5, 'LineStyle', '-'); % Box
                set(h(6,:), 'LineWidth', 2.5, 'LineStyle', '-'); % Median

                hold on;

                % Overlay individual data points (swarmchart)
                unique_groups = unique(grouping_data);
                for cl = unique_groups'
                    idx = (grouping_data == cl);
                    swarmchart(ones(sum(idx),1)*cl, plot_data(idx), 10, ...
                        'filled', 'MarkerFaceColor', plot_colors(cl,:), ...
                        'MarkerFaceAlpha', 0.7, 'XJitter', 'density', 'XJitterWidth', 0.4);
                end

                % Set axis properties
                ylabel(metric_name, 'Interpreter', 'none');
                xlabel('Genotype');
                xticklabels({'WT','OxtKO'});
                box off; axis square;

                %
                % Step 1: Calculate bin edges **once**, based on group 1 data
                [counts1, bin_edges] = histcounts(plot_data(grouping_data==1), 'BinMethod', 'auto');

                % Step 2: Plot first histogram
                subplot(2,2,4)
                histogram(plot_data(grouping_data==1), bin_edges, ...
                    'EdgeColor', 'none', 'FaceColor', plot_colors(1,:), 'Normalization', 'probability');
                hold on;

                % Step 3: Plot second histogram **using the same bin edges**
                histogram(plot_data(grouping_data==2), bin_edges, ...
                    'EdgeColor', 'none', 'FaceColor', plot_colors(2,:), 'Normalization', 'probability');
                hold on;

                % Step 4: Add median lines
                yl = ylim; % fix ylim before plotting lines
                if strcmp(myMethod,'median')
                    ll1 = line([median(plot_data(grouping_data==1)), median(plot_data(grouping_data==1))], yl);
                else
                    ll1 = line([mean(plot_data(grouping_data==1)), mean(plot_data(grouping_data==1))], yl);
                end
                ll1.Color = plot_colors(1,:);
                ll1.LineWidth = 2;

                if strcmp(myMethod,'median')
                    ll2 = line([median(plot_data(grouping_data==2)), median(plot_data(grouping_data==2))], yl);
                else
                    ll2 = line([mean(plot_data(grouping_data==2)), mean(plot_data(grouping_data==2))], yl);
                end
                ll2.Color = plot_colors(2,:);
                ll2.LineWidth = 2;

                % Final plot formatting
                xlabel(metric_name, 'Interpreter', 'none');
                ylabel('Probability');
                box off; axis square;

                % Add title with p-values
                sp = sgtitle(plot_title);
                sp.Interpreter = 'none';

                % Save Plot
                [annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
                exportgraphics(fig, fullfile(output_dir, sprintf('notBoxPlot_%s_by_Genotype_%s_%s.pdf', metric_name, day_range, myMethod)));

                % Save Source Data
                SourceDataTable = table(animal_ID,group, yData, 'VariableNames', {'Animal ID','Genotype','MetricValue'});
                writetable(SourceDataTable, fullfile(output_dir, sprintf('SourceData_%s_by_Genotype_%s_%s.csv', metric_name, day_range, myMethod)));

                % Save results for this metric
                if strcmp(myMethod,'mean')
                    Results = [Results; table({metric_name}, ...
                        AnalysisResults.diff_obs, AnalysisResults.pval_permutation, ...
                        AnalysisResults.pval_lme, AnalysisResults.pval_lme_log, ...
                        AnalysisResults.cohen_d, ...
                        AnalysisResults.skewness_WT, AnalysisResults.skewness_KO, AnalysisResults.skewness_perm, AnalysisResults.skewness_log, ...
                        'VariableNames', {'Metric','Observed_Diff','Pval_Perm','Pval_lme','Pval_lme_log','Cohens_d','Skewness_WT','Skewness_KO','Skewness_perm','Skewness_log'})];
                else
                    Results = [Results; table({metric_name}, ...
                        AnalysisResults.diff_obs, AnalysisResults.pval_permutation, ...
                        AnalysisResults.pval_lme, AnalysisResults.pval_lme_log, ...
                        AnalysisResults.cliffs_delta, ...
                        AnalysisResults.skewness_WT, AnalysisResults.skewness_KO, AnalysisResults.skewness_perm, AnalysisResults.skewness_log, ...
                        'VariableNames', {'Metric','Observed_Diff','Pval_Perm','Pval_lme','Pval_lme_log','Cliffs_d','Skewness_WT','Skewness_KO','Skewness_perm','Skewness_log'})];
                end
            end
            close all

            % Save results for this day range
            writetable(Results, fullfile(output_dir, sprintf('Genotype_vs_NoSeMazeMetrics_results_%s_%s.csv', day_range, myMethod)));
            disp(['Finished analysis for ', data_type, ' at ', day_range]);
        end
    end
end
disp('All genotype comparisons complete.');
