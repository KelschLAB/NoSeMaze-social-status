%% analyze_association_between_tube_metrics_and_manual_TT.m
% Author: Jonathan Reinwald
% Date: January 2025

% Purpose:
% Analyze associations between basic animal information (e.g., weight, age, genotype) 
% and NoSeMaze metrics (e.g., DSz, Rank, Wins/Losses). Includes scatter plots
% for each association and linear mixed-effects models accounting for Mouse_RFID.

%% Pre-Clearing
clear; clc; close all;

%% Set Directories
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory
processed_dir = fullfile(main_dir, 'data', 'processed');
results_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort', 'associations_between_metrics','tube','manualTT');
if ~isfolder(results_dir)
    mkdir(results_dir)
end
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');

% Add required paths
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));
addpath(genpath(fullfile(main_dir, 'src', 'analysis', 'associations_between_metrics')));

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
X_metrics.tube = {'DSz_Competition'};%, 'Rank_Competition', 'N_Wins_Competition', 'N_Losses_Competition', 'Fraction_Wins_Competition', 'Fraction_Losses_Competition', 'ELO_rand_Competition_Koptimal', 'ELO_nonrand_Competition_Koptimal',...
    % 'DSz_Chasing', 'Rank_Chasing', 'N_Wins_Chasing', 'N_Losses_Chasing', 'Fraction_Wins_Chasing', 'Fraction_Losses_Chasing', 'ELO_rand_Chasing_Koptimal', 'ELO_nonrand_Chasing_Koptimal'};
Y_metrics.tube = {'DSz_MTT'};%, , 'Rank_MTT', 'N_Wins_MTT', 'N_Losses_MTT', 'Fraction_Wins_MTT', 'Fraction_Losses_MTT',};
day_ranges = {'D1_21'}; % Select day_ranges of interest
LME_covariates = {'genotype'}; % Optional covariates: If none: LME_covariates = {''}

%% Iterate Through Data Types, NoSeMaze Metrics and Basic Metrics
for type_idx = 1:length(data_types)
    data_type = data_types{type_idx};
    
    % Filter input data (cohort selection)
    fields = fieldnames(summary_data.(data_type));
    for field_idx = 1:numel(fields)
        summary_data.(data_type).(fields{field_idx}) = summary_data.(data_type).(fields{field_idx})(ismember(summary_data.(data_type).(fields{field_idx}).cohort,cohortsTbl.cohort),:);
    end
    
    % Perform correlation matrix analysis
    correlation_output_dir = fullfile(results_dir, data_type, 'overview_association_matrix');
    if ~isfolder(correlation_output_dir)
        mkdir(correlation_output_dir);
    end
    
    % Correlation matrix is a good OVERVIEW
    generate_and_plot_association_matrix(summary_data.(data_type), X_metrics.(data_type), Y_metrics.(data_type), day_ranges, correlation_output_dir, 'Pearson', LME_covariates);

    for Y_metric_idx = 1:length(Y_metrics.(data_type))
        Y_metric = Y_metrics.(data_type){Y_metric_idx};

        % Create output directory
        Y_metric_output_dir = fullfile(results_dir, data_type, Y_metric);
        if ~isfolder(Y_metric_output_dir)
            mkdir(Y_metric_output_dir);
        end
        
        % Perform association analysis
        analyze_metric_associations(summary_data.(data_type), X_metrics.(data_type), Y_metric, day_ranges, Y_metric_output_dir, LME_covariates);
    end
end

disp('Analysis complete. Results saved.');
