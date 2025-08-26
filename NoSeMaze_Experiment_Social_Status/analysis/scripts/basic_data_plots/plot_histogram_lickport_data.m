%% plot_histogram_lickport_data.m
% Author: Jonathan Reinwald
% Date: January 2025
% Purpose: Generate and save plots histograms for pre-selected lickport data

%% Pre-Clearing
clear; clc; close all;

%% Set Directories
% Define main and subdirectories for data and results
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory
processed_dir = fullfile(main_dir, 'data', 'processed');
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');

% Add helper and plotting functions to the MATLAB path
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));
addpath(genpath(fullfile(main_dir, 'src', 'analysis', 'basic_data_plots')));

%% Load Summary Data
% Ensure the summary data file exists; load it if available
if ~isfile(summary_file)
    error('Summary file not found. Please run `main_build_cross_cohort_summary` first.');
end
load(summary_file, 'summary_data');

%% Cohort Selection
% Load cohort metadata and optionally filter cohorts
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));
% Example filter (uncomment if needed):
cohortsTbl = cohortsTbl(cohortsTbl.use_lickport == 1, :);

%% Define Parameters
% Specify data types, metrics, and day ranges to analyze
data_types = {'lickport'}; % Only analyze lickport metrics
metrics.lickport = {'baseline_rate_mean_omitfirst','cs_plus_modulation_peak_minus_base','cs_plus_modulation_peak',...
    'cs_minus_modulation_min_minus_base','cs_minus_modulation_min',...
    'correct_hit_rate','correct_hit_rate_allTrials',...
    'correct_rejection_rate','correct_rejection_rate_allTrials',...
    'delay_avoidance_learner',...
    'prediction_normalized_lick_diff','prediction_effectsize',...
    'cs_plus_ramping',...
    'cs_plus_detection_speed','cs_minus_detection_speed',...
    'cs_plus_valuation_peak','cs_plus_valuation_time_to_peak','cs_plus_valuation_peak_minus_base'...
    'baseline_rate_shaping','cs_plus_switch_latency_at_cs_shaping','cs_minus_switch_latency_at_cs_shaping','delay_avoidance_shaping',...
    };
day_ranges = {'D1_End'};

%% Main Loop: Iterate Through Data Types, Metrics, and Day Ranges
for type_idx = 1:length(data_types)
    data_type = data_types{type_idx};

    % Filter summary data to include only selected cohorts
    fields = fieldnames(summary_data.(data_type));
    for field_idx = 1:numel(fields)
        summary_data.(data_type).(fields{field_idx}) = ...
            summary_data.(data_type).(fields{field_idx})(...
            ismember(summary_data.(data_type).(fields{field_idx}).cohort, cohortsTbl.cohort), :);
    end

    % Iterate through day ranges
    for day_idx = 1:numel(day_ranges)
        range_name = day_ranges{day_idx};

        % Iterate through metrics
        for metric_idx = 1:numel(metrics.(data_type))
            current_metric = metrics.(data_type){metric_idx};

            % Create output directory for the cohort
            output_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort', 'basic_data_plots', data_type);
            if ~isfolder(output_dir)
                mkdir(output_dir);
            end

            % Data
            data = summary_data.(data_type).(range_name).(current_metric);

            % Calculate bin width dynamically using Freedman-Diaconis rule
            % Bin width formula: 2 * IQR(data) / (n^(1/3))
            IQR_data = iqr(data);            % Interquartile range
            n = length(data);                % Number of data points
            binWidth = 2 * IQR_data / (n^(1/3)); % Freedman-Diaconis rule for bin width

            % Round the bin width to make it more interpretable
            if binWidth < 0.5
                binWidth = 0.05; % Round to two decimal places for better clarity
            elseif binWidth > 0.5 && binWidth < 10
                binWidth = 0.5;
            elseif binWidth > 5 && binWidth < 50
                binWidth = 5;
            else
                disp('BinWidth not adapted');
            end

            % Calculate the number of bins based on the data range and bin width
            xMin = floor(min(data));
            xMax = ceil(max(data));
            binEdges = xMin:binWidth:xMax;

            % Plot histogram
            f = figure('Visible','on');
            H = histogram(data, binEdges, 'EdgeColor','none');

            % Set x-axis limits to be a little larger than the data range
            xlim([xMin, xMax]);

            % Title and labels
            tt=title(current_metric); tt.Interpreter = 'none';
            xlabel('Value');
            ylabel('Frequency');
                        
            % Save sorted data and ranks as a CSV
            DataTable = table(H.BinEdges(1:end-1)', H.BinEdges(2:end)', H.BinCounts(1:end)', ...
                'VariableNames', {'Edge1', 'Edge2', 'Count'});
            writetable(DataTable, fullfile(output_dir, ...
                ['SourceData_Histogram_' current_metric '_' range_name '.csv']));

            % Save the hierarchy score plot as a PDF
            exportgraphics(f, fullfile(output_dir, ...
                ['Histogram_' current_metric '_' range_name '.pdf']), ...
                'ContentType', 'vector', 'BackgroundColor', 'none');

            % Close all open figures to save memory
            % close all;

        end
    end
end





