%% main_build_cross_cohort_summary_lickport_data.m
% Author: Jonathan Reinwald
% Date: January 2025

% Purpose:
% Main script to build and update a cross-cohort summary structure (`summary_data`)
% for tube competitions and tube chasings, combining both into a single table per day range.

%% Pre-Clearing
clear all; clc;

%% Set Directories and Paths
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory
data_dir = fullfile(main_dir, 'data');
addpath(genpath(fullfile(main_dir, 'src', 'lickport', 'cross_cohort_summary')));
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));

% Path to cohort configuration files
config_dir = fullfile(main_dir, 'config');

% Path to processed data
processed_dir = fullfile(data_dir, 'processed');
summary_file_dir = fullfile(processed_dir, 'cross_cohort_files');
if ~isfolder(summary_file_dir)
    mkdir(summary_file_dir);
end

% Path to save the summary data
summary_file = fullfile(summary_file_dir, 'summary_data.mat');
data_transformation_file = fullfile(summary_file_dir, 'data_transformation.mat');

%% Load Cohort Information
cohortsTbl = readtable(fullfile(config_dir, 'cohorts_info.csv')); % Load cohort metadata
cohortsDetails = readtable(fullfile(config_dir, 'cohorts_details.csv')); % Load animal details

%% Initialize or Load Existing Summary Data
if isfile(summary_file)
    load(summary_file, 'summary_data'); % Load existing summary structure
else
    summary_data = struct(); % Initialize empty structure
end

%% Define Day Ranges
day_ranges = {'D1_End'};

%% Load cross_cohort_files
data_dir = fullfile(processed_dir, 'cross_cohort_files','lickport');
lick_params_file = fullfile(data_dir, 'lick_params_cross_cohort.mat');

% Skip if the file does not exist
if ~isfile(lick_params_file)
    error('No cross-cohort lick params file not found.');
end

% Load the hierarchy data
load(lick_params_file, 'lick_params');

% Filter animal details for the current cohort
cohortDetails = cohortsDetails;

% Create lickport table
lickport_table = struct2table(lick_params);

% Change ID name to Mouse_RFID
lickport_table.Properties.VariableNames{find(strcmp(lickport_table.Properties.VariableNames,'ID'))} = 'Mouse_RFID';

% Selection of metrics
plot_features = {'baseline_rate_mean_omitfirst','baseline_last_trial_rewarded_omitfirst','baseline_last_trial_non_rewarded_omitfirst',... % pre-CS lick rates
    'impulsivity_trialStart','impulsivity_odorStart','impulsivity_rewardStart',...
    'impulsivity_trialStart_norm','impulsivity_odorStart_norm','impulsivity_rewardStart_norm',...
    'prediction_normalized_lick_diff','prediction_lick_diff','prediction_effectsize','prediction_significance',... % prediction of reward based on previous trial's outcome
    'baseline_rate_shaping','cs_plus_switch_latency_at_cs_shaping','cs_plus_switch_latency_at_us_shaping',... % shaping
    'cs_minus_switch_latency_at_cs_shaping','cs_minus_switch_latency_at_us_shaping','delay_avoidance_shaping',... % shaping
    'pause_duration_at_CS_shaping_rev_1toLast','pause_duration_at_US_shaping_rev_1toLast',... % shaping
    'cs_plus_modulation_peak_minus_base','cs_plus_modulation_peak','cs_plus_modulation_averaged',...% CS modulation
    'cs_minus_modulation_min_minus_base','cs_minus_modulation_min','cs_minus_modulation_averaged',...
    'cs_plus_ramping',...
    'cs_plus_detection_speed','cs_minus_detection_speed',...
    'cs_plus_valuation_peak','cs_plus_valuation_peak_minus_base','cs_plus_valuation_time_to_peak',...
    'cs_plus_switch_latency_at_cs_rev1','cs_plus_switch_latency_at_cs_rev2','cs_plus_switch_latency_at_cs_rev3',... % switch latencies
    'cs_plus_switch_latency_at_cs_rev4','cs_plus_switch_latency_at_cs_mean','cs_plus_switch_latency_at_cs_median',...
    'cs_minus_switch_latency_at_cs_rev1','cs_minus_switch_latency_at_cs_rev2','cs_minus_switch_latency_at_cs_rev3',...
    'cs_minus_switch_latency_at_cs_rev4','cs_minus_switch_latency_at_cs_mean','cs_minus_switch_latency_at_cs_median',...
    'cs_plus_switch_latency_at_us_rev1','cs_plus_switch_latency_at_us_rev2','cs_plus_switch_latency_at_us_rev3',...
    'cs_plus_switch_latency_at_us_rev4','cs_plus_switch_latency_at_us_mean','cs_plus_switch_latency_at_us_median',...
    'cs_minus_switch_latency_at_us_rev1','cs_minus_switch_latency_at_us_rev2','cs_minus_switch_latency_at_us_rev3',...
    'cs_minus_switch_latency_at_us_rev4','cs_minus_switch_latency_at_us_mean','cs_minus_switch_latency_at_us_median',...
    'pause_duration_at_CS_rev1','pause_duration_at_CS_rev2','pause_duration_at_CS_rev3','pause_duration_at_CS_rev4',... % pauses
    'pause_duration_at_US_rev1','pause_duration_at_US_rev2','pause_duration_at_US_rev3', 'pause_duration_at_US_rev4',...
    'delay_avoidance_learner',...% delay avoidance
    'giving_up_at_CS_rev1','giving_up_at_CS_rev2','giving_up_at_CS_rev3','giving_up_at_CS_rev4',... % giving up
    'giving_up_at_US_rev1','giving_up_at_US_rev2','giving_up_at_US_rev3','giving_up_at_US_rev4',...
    'correct_rejection_rate','correct_rejection_rate_allTrials',... % correct rejection rates
    'correct_rejection_rate_rev1','correct_rejection_rate_rev2','correct_rejection_rate_rev3','correct_rejection_rate_rev4',...
    'correct_hit_rate','correct_hit_rate_allTrials',... % correct hit rates
    'correct_hit_rate_rev1','correct_hit_rate_rev2','correct_hit_rate_rev3','correct_hit_rate_rev4',...
    'cs_minus_time_to_criterion_phase1','cs_minus_time_to_criterion_phase2','cs_minus_time_to_criterion_phase3','cs_minus_time_to_criterion_phase4','cs_minus_time_to_criterion_phase5','cs_minus_time_to_criterion_phase6','cs_minus_time_to_criterion_mean','cs_minus_time_to_criterion_median',...
    'cs_plus_time_to_criterion_phase1','cs_plus_time_to_criterion_phase2','cs_plus_time_to_criterion_phase3','cs_plus_time_to_criterion_phase4','cs_plus_time_to_criterion_phase5','cs_plus_time_to_criterion_phase6','cs_plus_time_to_criterion_mean','cs_plus_time_to_criterion_median'}

myLearningTable=lickport_table(:,[{'Mouse_RFID'},{'cohort'},plot_features(:)']);

plot_output = fullfile(main_dir,'/results/figures/cross_cohort/basic_data_checks/lickport/histogram_normal_distribution');
if ~isfolder(plot_output); mkdir(plot_output); end

% plot distributions of selected features
if 1==0
    for ii = 1:numel(plot_features)
        f=figure('Visible','off');
        clear plot_data
        plot_data = lickport_table.(plot_features{ii});
        histogram(plot_data(~isnan(plot_data)),'EdgeColor','none');%,'BinWidth',BinWidth_features{ii});
        box off
        ax=gca;
        % ax.XLim=XLimits_features{ii};
        title(plot_features{ii},'Interpreter','none');
        saveas(f,fullfile(plot_output,[plot_features{ii},'_original.png']),'png');
        saveas(f,fullfile(plot_output,[plot_features{ii},'_original.pdf']),'pdf');
        close all;
    end
end

% boxcox transformation and plotting
for ii = 1:numel(plot_features)

    % preparation
    clear myData
    myData = myLearningTable.(plot_features{ii});
    myData_red = myData(~isnan(myData) & ~isinf(myData));
    data_transformation(ii).metric_name = plot_features{ii};
    data_transformation(ii).original_data = myData_red;

    % 1. test for normality with sw and ks
    [data_transformation(ii).H_sw, data_transformation(ii).pValue_sw, data_transformation(ii).W_sw] = swtest(myData_red);
    [data_transformation(ii).H_ks, data_transformation(ii).pValue_ks, ~, ~] = kstest(myData_red);

    % 2. data transformation
    if data_transformation(ii).H_sw==1 || data_transformation(ii).H_ks==1

        % 2.1 tranformation (searching for optimal lambda and transform
        % data)
        % Add a small constant to avoid zero or negative values
        if min(myData_red)<0
            epsilon = ceil(abs(min(myData_red)));
        else
            epsilon = 1e-3; % Small constant
        end
        data_adj = myData_red + epsilon;

        % Define a range of lambda values to test
        lambda_values = -10:0.01:10;
        % Initialize variable to store the maximum log-likelihood and the best lambda
        max_log_likelihood = -Inf;
        best_lambda = NaN;
        % Iterate over each lambda value
        for lambda = lambda_values
            if lambda == 0
                % For lambda = 0, use log transformation
                transformed_data = log(data_adj);
            else
                % Box-Cox transformation for lambda != 0
                transformed_data = (data_adj .^ lambda - 1) / lambda;
            end
            % Calculate the log-likelihood
            n = length(data_adj);
            log_likelihood = -n/2 * log(var(transformed_data)) + (lambda - 1) * sum(log(data_adj));
            % Update the best lambda if the current log-likelihood is higher
            if log_likelihood > max_log_likelihood
                max_log_likelihood = log_likelihood;
                best_lambda = lambda;
            end
        end
        % Display the optimal lambda
        disp('Optimal Lambda:');
        disp(best_lambda);

        if best_lambda == 0
            transformed_data_boxcox = myData_red;
            transformed_data_log10 = log10(data_adj);
            transformed_data_cuberoot = data_adj .^ (1/3);
        else
            transformed_data_boxcox = (data_adj .^ best_lambda - 1) / best_lambda;
            transformed_data_log10 = log10(data_adj);
            transformed_data_cuberoot = data_adj .^ (1/3);
        end

        % 2.2 test for normality with transformed data
        [data_transformation(ii).H_sw_transformed_boxcox, data_transformation(ii).pValue_sw_transformed_boxcox, ~] = swtest(transformed_data_boxcox);
        [data_transformation(ii).H_ks_transformed_boxcox, data_transformation(ii).pValue_ks_transformed_boxcox, ~, ~] = kstest(transformed_data_boxcox);

        [data_transformation(ii).H_sw_log10, data_transformation(ii).pValue_sw_log10, ~] = swtest(transformed_data_log10);
        [data_transformation(ii).H_ks_log10, data_transformation(ii).pValue_ks_log10, ~, ~] = kstest(transformed_data_log10);

        [data_transformation(ii).H_sw_cuberoot, data_transformation(ii).pValue_sw_cuberoot, ~] = swtest(transformed_data_cuberoot);
        [data_transformation(ii).H_ks_cuberoot, data_transformation(ii).pValue_ks_cuberoot, ~, ~] = kstest(transformed_data_cuberoot);

        [data_transformation(ii).skewness_original] = skewness(data_adj);
        [data_transformation(ii).skewness_cuberoot] = skewness(transformed_data_cuberoot);
        [data_transformation(ii).skewness_log10] = skewness(transformed_data_log10);
        [data_transformation(ii).skewness_boxcox] = skewness(transformed_data_boxcox);

        data_transformation(ii).transformed_data_boxcox = transformed_data_boxcox;
        data_transformation(ii).transformed_data_log10 = transformed_data_log10;
        data_transformation(ii).transformed_data_cuberoot = transformed_data_cuberoot;
        data_transformation(ii).best_lambda = best_lambda;

        % 2.3 redo plotting
        f=figure('Visible','off');
        histogram(data_transformation(ii).transformed_data_boxcox,'EdgeColor','none');
        box off
        title([plot_features{ii} ' (transformed boxcox)'],'Interpreter','none');
        saveas(f,fullfile(plot_output,[plot_features{ii},'_transformedBoxCox.png']),'png');
        saveas(f,fullfile(plot_output,[plot_features{ii},'_transformedBoxCox.pdf']),'pdf');
        close all;

        f=figure('Visible','off');
        histogram(data_transformation(ii).transformed_data_log10,'EdgeColor','none');
        box off
        title([plot_features{ii} ' (transformed log10)'],'Interpreter','none');
        saveas(f,fullfile(plot_output,[plot_features{ii},'_transformedLog10.png']),'png');
        saveas(f,fullfile(plot_output,[plot_features{ii},'_transformedLog10.pdf']),'pdf');
        close all;

        f=figure('Visible','off');
        histogram(data_transformation(ii).transformed_data_cuberoot,'EdgeColor','none');
        box off
        title([plot_features{ii} ' (transformed cuberoot)'],'Interpreter','none');
        saveas(f,fullfile(plot_output,[plot_features{ii},'_transformedcuberoot.png']),'png');
        saveas(f,fullfile(plot_output,[plot_features{ii},'_transformedcuberoot.pdf']),'pdf');
        close all;
    end

    % preparation
    myLearningTable = [myLearningTable,table(nan(height(myLearningTable),1),'VariableNames',{['log10_' plot_features{ii}]})];
    myLearningTable.(['log10_' plot_features{ii}])(find(~isnan(myData) & ~isinf(myData))) = data_transformation(ii).transformed_data_log10;

    myLearningTable = [myLearningTable,table(nan(height(myLearningTable),1),'VariableNames',{['boxcox_' plot_features{ii}]})];
    myLearningTable.(['boxcox_' plot_features{ii}])(find(~isnan(myData) & ~isinf(myData))) = data_transformation(ii).transformed_data_boxcox;

    myLearningTable = [myLearningTable,table(nan(height(myLearningTable),1),'VariableNames',{['cuberoot_' plot_features{ii}]})];
    myLearningTable.(['cuberoot_' plot_features{ii}])(find(~isnan(myData) & ~isinf(myData))) = data_transformation(ii).transformed_data_cuberoot ;
end

% [r,p]=corr(table2array(myLearningTable(:,3:77)),'rows','pairwise','type','Spearman'); names=myLearningTable.Properties.VariableNames(3:77);
% figure; imagesc(r.*(p<0.05)); ax=gca;
% ax.XTick = 1:length(names);
% ax.YTick = 1:length(names);
% ax.YTickLabel = names; ax.XTickLabel = names;
% ax.TickLabelInterpreter = 'none';
% ax.CLim=[-1,1]
% ax.Colormap = jet

% Align metrics with cohortDetails based on Mouse_RFID
full_table = outerjoin(cohortDetails, myLearningTable, ...
    'Keys', {'Mouse_RFID','cohort'}, 'MergeKeys', true, 'Type', 'Left');

summary_data.lickport.(day_ranges{1}) = full_table;

%% Save Updated Summary Data
save(summary_file, 'summary_data');
fprintf('Summary data saved to %s/n', summary_file);
save(data_transformation_file, 'data_transformation');
fprintf('Data transformation information saved to %s/n', data_transformation_file);

figure; histogram([data_transformation.skewness_original],'BinWidth',0.5,'EdgeColor','none');
hold on; histogram([data_transformation.skewness_boxcox],'BinWidth',0.5,'EdgeColor','none');
hold on; histogram([data_transformation.skewness_log10],'BinWidth',0.5,'EdgeColor','none');
hold on; histogram([data_transformation.skewness_cuberoot],'BinWidth',0.5,'EdgeColor','none');

myFieldNames = fieldnames(data_transformation);
myFieldNames = myFieldNames(contains(myFieldNames,'skewness'));
p_table = table();
t_table = table();

for ix = 1:length(myFieldNames)
    for jx = ix+1:length(myFieldNames)
        [clusters, p_values, t_sums, permutation_distribution ] = permutest( [data_transformation.(myFieldNames{ix})],[data_transformation.(myFieldNames{jx})],true, ...
            0.05, 10000, true);
        p_table.([myFieldNames{ix} '_vs_' myFieldNames{jx} '_p']) = p_values;
        if isempty(t_sums)
            t_table.([myFieldNames{ix} '_vs_' myFieldNames{jx} '_tsums']) = nan;
        else
            t_table.([myFieldNames{ix} '_vs_' myFieldNames{jx} '_tsums']) = t_sums;
        end
        disp(['p: ' num2str(p_values) '; t: ' num2str(t_sums)])
    end
end
