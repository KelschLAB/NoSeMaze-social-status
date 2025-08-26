%% analyze_association_between_modalities.m
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
X_metrics.tube =  {'DSz_Competition', 'N_Events_Competition',...
    'N_Wins_Competition','Fraction_Wins_Competition', 'N_Losses_Competition', 'Fraction_Losses_Competition',...
    'DSz_Chasing', 'N_Events_Chasing',...
    'N_Wins_Chasing','Fraction_Wins_Chasing', 'N_Losses_Chasing', 'Fraction_Losses_Chasing'};

% {'DSz_Competition', 'Fraction_Wins_Chasing', 'Fraction_Losses_Chasing'};
Y_metrics.tube =  {'DSz_Competition', 'N_Events_Competition',...
    'N_Wins_Competition','Fraction_Wins_Competition', 'N_Losses_Competition', 'Fraction_Losses_Competition',...
    'DSz_Chasing', 'N_Events_Chasing',...
    'N_Wins_Chasing','Fraction_Wins_Chasing', 'N_Losses_Chasing', 'Fraction_Losses_Chasing'};

% {'DSz_Competition', 'Fraction_Wins_Chasing', 'Fraction_Losses_Chasing'};
X_metrics.lickport = {'correct_hit_rate_allTrials','correct_rejection_rate_allTrials',...
    'baseline_rate_mean_omitfirst','cs_plus_modulation_peak',...
    'cs_plus_switch_latency_at_cs_rev1','cs_minus_switch_latency_at_cs_rev1','cs_plus_switch_latency_at_cs_rev2','cs_minus_switch_latency_at_cs_rev2',...
    'cs_plus_switch_latency_at_cs_rev3','cs_minus_switch_latency_at_cs_rev3','cs_plus_switch_latency_at_cs_rev4','cs_minus_switch_latency_at_cs_rev4','cs_plus_switch_latency_at_cs_median','cs_minus_switch_latency_at_cs_median'};
Y_metrics.lickport = {'correct_hit_rate_allTrials','correct_rejection_rate_allTrials',...
    'baseline_rate_mean_omitfirst','cs_plus_modulation_peak',...
    'cs_plus_switch_latency_at_cs_rev1','cs_minus_switch_latency_at_cs_rev1','cs_plus_switch_latency_at_cs_rev2','cs_minus_switch_latency_at_cs_rev2',...
    'cs_plus_switch_latency_at_cs_rev3','cs_minus_switch_latency_at_cs_rev3','cs_plus_switch_latency_at_cs_rev4','cs_minus_switch_latency_at_cs_rev4','cs_plus_switch_latency_at_cs_median','cs_minus_switch_latency_at_cs_median'};
% X_metrics.lickport = {'cuberoot_correct_hit_rate_allTrials','cuberoot_correct_rejection_rate_allTrials',...
%     'cuberoot_baseline_rate_mean_omitfirst','cuberoot_cs_plus_modulation_peak','cuberoot_cs_minus_modulation_min',...
%     'cuberoot_cs_plus_switch_latency_at_cs_rev1','cuberoot_cs_minus_switch_latency_at_cs_rev1','cuberoot_cs_plus_switch_latency_at_cs_rev2','cuberoot_cs_minus_switch_latency_at_cs_rev2',...
%     'cuberoot_cs_plus_switch_latency_at_cs_rev3','cuberoot_cs_minus_switch_latency_at_cs_rev3','cuberoot_cs_plus_switch_latency_at_cs_rev4','cuberoot_cs_minus_switch_latency_at_cs_rev4'};
% Y_metrics.lickport = {'cuberoot_correct_hit_rate_allTrials','cuberoot_correct_rejection_rate_allTrials',...
%     'cuberoot_baseline_rate_mean_omitfirst','cuberoot_cs_plus_modulation_peak','cuberoot_cs_minus_modulation_min',...
%     'cuberoot_cs_plus_switch_latency_at_cs_rev1','cuberoot_cs_minus_switch_latency_at_cs_rev1','cuberoot_cs_plus_switch_latency_at_cs_rev2','cuberoot_cs_minus_switch_latency_at_cs_rev2',...
%     'cuberoot_cs_plus_switch_latency_at_cs_rev3','cuberoot_cs_minus_switch_latency_at_cs_rev3','cuberoot_cs_plus_switch_latency_at_cs_rev4','cuberoot_cs_minus_switch_latency_at_cs_rev4'};


day_ranges = {'D1_End'}; % Select day_ranges of interest
LME_covariates = {'genotype'}; % Optional covariates: If none: LME_covariates = {''}

    
% Output directory
correlation_output_dir = fullfile(results_dir, 'overview_association_matrix');
if ~isfolder(correlation_output_dir)
    mkdir(correlation_output_dir);
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
    % Correlation matrix is a good OVERVIEW
    generate_and_plot_association_matrix(data_table_combined, X_metrics_combined.(range_name), Y_metrics_combined.(range_name), day_ranges, correlation_output_dir, 'Spearman');
    
end

% 
% 
% 
% 
%         for Y_metric_idx = 1:length(Y_metrics.(data_type))
%             Y_metric = Y_metrics.(data_type){Y_metric_idx};
% 
%             % Create output directory
%             Y_metric_output_dir = fullfile(results_dir, data_type, Y_metric);
%             if ~isfolder(Y_metric_output_dir)
%                 mkdir(Y_metric_output_dir);
%             end
% 
%             % Perform association analysis
%             analyze_metric_associations(summary_data.(data_type), X_metrics.(data_type), Y_metric, day_ranges.(data_type), Y_metric_output_dir, LME_covariates);
%         end
%     end
% end

disp('Analysis complete. Results saved.');
