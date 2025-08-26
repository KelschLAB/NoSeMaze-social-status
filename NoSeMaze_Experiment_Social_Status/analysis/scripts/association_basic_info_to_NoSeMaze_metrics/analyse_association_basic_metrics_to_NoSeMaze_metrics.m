%% analyze_association_basic_info_to_NoSeMaze_metrics.m
% Author: Jonathan Reinwald
% Date: January 2025

% Purpose:
% Analyze associations between basic animal information (e.g., weight, age, genotype) 
% and NoSeMaze metrics (e.g., DSz, Rank, Wins/Losses). Includes scatter plots
% for each association and linear mixed-effects models accounting for Mouse_RFID.

%% Pre-Clearing
clear; clc; close all;

%% Set Directories
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory
processed_dir = fullfile(main_dir, 'data', 'processed');
results_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort', 'associations_basic_info_to_NoSeMaze_metrics');
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');

% Add required paths
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));
addpath(genpath(fullfile(main_dir, 'src', 'analysis', 'associations_basic_info_to_NoSeMaze_metrics')));
% addpath(genpath(fullfile(main_dir, 'src', 'analysis', 'associations_between_metrics')));

%% Load Summary Data
if ~isfile(summary_file)
    error('Summary file not found. Please run `main_build_cross_cohort_summary` first.');
end
load(summary_file, 'summary_data');

%% Cohort Selection
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));
cohortsTbl = cohortsTbl(cohortsTbl.use_tube == 1, :);

%% Define Parameters
data_types = {'tube'}; % Analyze tube metrics
Basic_metrics = {'weight_pre_zscore'};%'weight_pre', 'weight_pre_zscore', 'weight_diff', 'age', 'age_zscore', 'genotype'}; % Binary & continuous
% NoSeMaze_metrics.tube = {'DSz_Competition','cuberoot_Fraction_Wins_Competition', 'cuberoot_Fraction_Losses_Competition','Fraction_Wins_Chasing', 'Fraction_Losses_Chasing', 'N_Events_Chasing'};
NoSeMaze_metrics.tube = {'cuberoot_Fraction_Wins_Chasing', 'cuberoot_Fraction_Losses_Chasing'};
% , 'N_Events_Competition', 'Fraction_Wins_Competition', 'Fraction_Losses_Competition', 'N_Wins_Competition', 'N_Losses_Competition', ...
%     'detection_count', 'detection_percentage', 'detections_per_day', 'DSz_Chasing', 'N_Events_Chasing', 'cuberoot_Fraction_Wins_Chasing', 'cuberoot_Fraction_Losses_Chasing'};
day_ranges.tube = {'D1_21'}; % Select day_ranges of interest
NoSeMaze_metrics.lickport = {'correct_hit_rate_allTrials','correct_rejection_rate_allTrials',...
    'baseline_rate_mean_omitfirst','cs_plus_modulation_peak',...
    'cs_plus_switch_latency_at_cs_rev1','cs_minus_switch_latency_at_cs_rev1','cs_plus_switch_latency_at_cs_rev2','cs_minus_switch_latency_at_cs_rev2',...
    'cs_plus_switch_latency_at_cs_rev3','cs_minus_switch_latency_at_cs_rev3','cs_plus_switch_latency_at_cs_rev4','cs_minus_switch_latency_at_cs_rev4','cs_plus_switch_latency_at_cs_median','cs_minus_switch_latency_at_cs_median'};
% {'baseline_rate_mean_omitfirst','cs_plus_modulation_peak_minus_base','cs_plus_modulation_peak',...
%     'cs_minus_modulation_min_minus_base','cs_minus_modulation_min',...
%     'correct_hit_rate','correct_hit_rate_allTrials',...
%     'correct_rejection_rate','correct_rejection_rate_allTrials',...
%     'delay_avoidance_learner',...
%     'prediction_normalized_lick_diff','prediction_effectsize',...
%     'cs_plus_ramping',...
%     'cs_plus_detection_speed','cs_minus_detection_speed',...
%     'cs_plus_valuation_peak','cs_plus_valuation_time_to_peak','cs_plus_valuation_peak_minus_base'...
%     'baseline_rate_shaping','cs_plus_switch_latency_at_cs_shaping','cs_minus_switch_latency_at_cs_shaping','delay_avoidance_shaping',...
%     };
day_ranges.lickport = {'D1_21'}; % Select day_ranges of interest
LME_covariates = {'genotype'}; % {'repetition'};

%% Iterate Through Data Types, NoSeMaze Metrics and Basic Metrics
for type_idx = 1%:length(data_types)
    data_type = data_types{type_idx};
    
    % Filter input data (cohort selection)
    fields = fieldnames(summary_data.(data_type));
    for field_idx = 1:numel(fields)
        summary_data.(data_type).(fields{field_idx}) = summary_data.(data_type).(fields{field_idx})(ismember(summary_data.(data_type).(fields{field_idx}).cohort,cohortsTbl.cohort),:);
    end
    
    %% Standardization of age and weight
    % Assuming your table is named 'summary_data.tube.D1_21'
    tbl = summary_data.(data_type).(day_ranges.(data_type){1});

    % Get unique cohorts
    cohorts = unique(tbl.cohort);

    % Initialize new variables
    tbl.weight_pre_zscore = NaN(height(tbl), 1);
    tbl.weight_post_zscore = NaN(height(tbl), 1);
    tbl.age_zscore = NaN(height(tbl), 1);
    tbl.weight_diff_percentual = tbl.weight_diff./tbl.weight_pre;

    % Loop through each cohort and standardize variables
    for i = 1:length(cohorts)
        cohort_idx = strcmp(tbl.cohort, cohorts{i}); % Get index for cohort
        cohort_weight_post = tbl.weight_post(cohort_idx);
        cohort_weight_pre = tbl.weight_pre(cohort_idx);
        cohort_age = tbl.age(cohort_idx);

        % Compute z-scores
        tbl.weight_pre_zscore(cohort_idx) = (cohort_weight_pre - mean(cohort_weight_pre)) ./ std(cohort_weight_pre);
        tbl.age_zscore(cohort_idx) = (cohort_age - mean(cohort_age)) ./ std(cohort_age);
        tbl.weight_post_zscore(cohort_idx) = (cohort_weight_post - mean(cohort_weight_post)) ./ std(cohort_weight_post);
    end

    % Assign back to the main table
    summary_data.(data_type).(day_ranges.(data_type){1}) = tbl;
    %%

    % Perform correlation matrix analysis
    correlation_output_dir = fullfile(results_dir, data_type, 'overview_association_matrix');
    if ~isfolder(correlation_output_dir)
        mkdir(correlation_output_dir);
    end
    
    % Association matrix is a good OVERVIEW
    generate_and_plot_association_matrix_basic_metrics(summary_data.(data_type), Basic_metrics, NoSeMaze_metrics.(data_type), day_ranges.(data_type), correlation_output_dir, 'Spearman');

    for NoSeMaze_metric_idx = 1:length(NoSeMaze_metrics.(data_type))
        NoSeMaze_metric = NoSeMaze_metrics.(data_type){NoSeMaze_metric_idx};

        % Create output directory
        NoSeMaze_metric_output_dir = fullfile(results_dir, data_type, NoSeMaze_metric);
        if ~isfolder(NoSeMaze_metric_output_dir)
            mkdir(NoSeMaze_metric_output_dir);
        end
        
        % Perform association analysis
        analyze_basic_metric_associations(summary_data.(data_type), NoSeMaze_metric, Basic_metrics, day_ranges.(data_type), NoSeMaze_metric_output_dir, LME_covariates);
    end
end

disp('Analysis complete. Results saved.');
