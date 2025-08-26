%% check_stability_over_time.m
% Author: Jonathan Reinwald
% Date: January 2025

% This script analyzes the stability of metrics (e.g., DSz, Rank, Wins/Losses)
% across different day ranges for various data types (e.g., tube, lickport).

%% Pre-Clearing
clear; clc; close all;

%% Set Directories
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory
processed_dir = fullfile(main_dir, 'data', 'processed');
results_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort', 'stability_over_time');
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');

% Add required paths
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));
addpath(genpath(fullfile(main_dir, 'src', 'analysis', 'stability_over_time')));

%% Load Summary Data
if ~isfile(summary_file)
    error('Summary file not found. Please run `main_build_cross_cohort_summary` first.');
end
load(summary_file, 'summary_data');

%% Cohort Selection
% Filter the table to include only cohorts with tube data of interest (use_tube == 1)
cohortsTbl = readtable(fullfile(main_dir,'config','cohorts_info.csv'));
cohortsTbl = cohortsTbl(cohortsTbl.use_tube==1,:);

%% Define Data Types, Metrics, and Day Ranges
data_types = {'tube'}; %{'tube', 'lickport'}; % Extendable for additional data types
% metrics.tube = {'boxcox_N_Wins_Chasing', 'boxcox_N_Losses_Chasing','boxcox_Fraction_Wins_Chasing', 'boxcox_Fraction_Losses_Chasing','boxcox_N_Events_Chasing'};
metrics.tube = {'DSz_Competition','cuberoot_Fraction_Wins_Chasing','cuberoot_Fraction_Losses_Chasing','cuberoot_N_Events_Chasing',...
    'Fraction_Wins_Chasing','Fraction_Losses_Chasing','N_Events_Chasing'};
% metrics.tube = {'Fraction_Wins_Chasing', 'Fraction_Losses_Chasing','N_Events_Chasing','N_Wins_Chasing', 'N_Losses_Chasing'};
% 'DSz_Competition', 'Rank_Competition', 'N_Wins_Competition', ...
%                 'N_Losses_Competition', 'Fraction_Wins_Competition', 'Fraction_Losses_Competition', ...
%                 'DSz_Chasing', 'Rank_Chasing', 
%metrics.lickport = {'test'}; % Example lickport metrics
day_ranges = {'D1_7','D8_14','D15_21'}; % Select day_ranges of interest
LME_covariates = {'genotype'};%{'age','genotype'}; % Option to integrate covariates into the LME; note that the RFID is anyway considered as a random effects grouping variable

%% Generate Figures and Results for Each Data Type
for type_idx = 1:length(data_types)
    data_type = data_types{type_idx};
    
    for metric_idx = 1:length(metrics.(data_type))
        metric = metrics.(data_type){metric_idx}; % e.g., 'DSz_Competition'
        
        % Create output directory for the current metric
        metric_output_dir = fullfile(results_dir, data_type, metric);
        if ~isfolder(metric_output_dir)
            mkdir(metric_output_dir);
        end
        
        % Filter input data
        fields = fieldnames(summary_data.(data_type));
        for field_idx = 1:numel(fields)
            summary_data.(data_type).(fields{field_idx}) = summary_data.(data_type).(fields{field_idx})(ismember(summary_data.(data_type).(fields{field_idx}).cohort,cohortsTbl.cohort),:);
        end
        
        % Generate plots and perform LME analysis
        plot_and_analyze_metric(summary_data.(data_type), metric, day_ranges, metric_output_dir,LME_covariates);
    end
end

disp('Analysis complete. Results saved.');

