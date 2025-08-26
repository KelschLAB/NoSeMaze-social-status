%% main_build_cross_cohort_summary_tube_data.m
% Author: Jonathan Reinwald
% Date: January 2025

% Purpose:
% Main script to build and update a cross-cohort summary structure (`summary_data`)
% for tube competitions and tube chasings, combining both into a single table per day range.

%% Pre-Clearing
clear; clc;

%% Set Directories and Paths
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory
data_dir = fullfile(main_dir, 'data');
addpath(genpath(fullfile(main_dir, 'src', 'tube', 'cross_cohort_summary')));
addpath(genpath(fullfile(main_dir, 'src', 'tube', 'preprocessing')));
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
day_ranges = {'D1_21', 'D1_7', 'D8_14', 'D15_21', 'D1_14', 'D1_End', 'Last14', };

%% Process Data of all Cohorts across all Day Ranges

% Iterate over day ranges
for range_idx = 1:length(day_ranges)
    range_name = day_ranges{range_idx}; % e.g., 'D1_7'

    % Initialize a combined table for competitions and chasings (Outside of cohort loop)
    tube_table = [];

    % Process each cohort
    for cohort_idx = 1:height(cohortsTbl)       
        cohort = cohortsTbl.cohort{cohort_idx};
        disp(['Range: ' range_name ', Processing ' cohort ]);
        cohort_data_dir = fullfile(processed_dir, cohort, 'tube', 'full_hierarchy_files');
        full_hierarchy_file = fullfile(cohort_data_dir, ['full_hierarchy_' cohort '.mat']);

        % Skip if the file does not exist
        if ~isfile(full_hierarchy_file)
            fprintf('Skipping cohort %s: hierarchy file not found./n', cohort);
            continue;
        end

        % Load the hierarchy data
        load(full_hierarchy_file, 'full_hierarchy');

        % Determine the number of days for this cohort
        num_days = length(full_hierarchy);  % Adjust according to your data structure

        % Define day range (start and end days based on the current day range)
        if contains(range_name, 'End')
            % For ranges with "End"
            range_parts = regexp(range_name, 'D(\d+)_', 'tokens'); % Extract start day
            start_day = str2double(range_parts{1}{1});
            end_day = num_days;
        elseif strcmp(range_name, 'Last14')
            % 'Last14' is the last 14 days (or fewer if there are fewer than 14 days)
            start_day = max(num_days - 13, 1); % Ensure we don't go before day 1
            end_day = num_days; % Go to the last available day
        else
            % For other ranges like 'D1_7', 'D8_14', etc., extract the start and end from the string
            range_parts = regexp(range_name, 'D(\d+)_(\d+)', 'tokens'); % Extract start and end from 'D1_7'

            if ~isempty(range_parts)
                start_day = str2double(range_parts{1}{1}); % Convert to number
                end_day = str2double(range_parts{1}{2}); % Convert to number
                % Ensure the end day does not exceed the available days
                end_day = min(end_day, num_days);
            else
                continue; % Skip if no valid range found
            end
        end

        % Definition of day_range
        day_range = start_day:end_day;

        % Filter animal details for the current cohort
        cohortDetails = cohortsDetails(strcmp(cohortsDetails.cohort, cohort), :);

        % Plot output to check normality 
        plot_output = fullfile(main_dir,'results','figures','cross_cohort','basic_data_checks','tube',range_name,'histogram_normal_distribution');
        if ~isfolder(plot_output) 
            mkdir(plot_output); 
        end

        % Process tube competitions and append suffix '_Competition'
        competition_table = process_tube_data(...
            full_hierarchy, day_range, cohortDetails, 'match_matrix', '_Competition',plot_output);

        % Process tube chasings and append suffix '_Chasing'
        chasing_table = process_tube_data(...
            full_hierarchy, day_range, cohortDetails, 'match_matrix_chasing', '_Chasing',plot_output);

        % Dynamically identify shared fields (excluding the key 'Mouse_RFID')
        shared_fields = intersect(competition_table.Properties.VariableNames, ...
            chasing_table.Properties.VariableNames);
        shared_fields = setdiff(shared_fields, {'Mouse_RFID'}); % Exclude the key

        % Remove shared fields from chasing_table
        chasing_table = removevars(chasing_table, shared_fields);

        % Merge the tables using outerjoin
        tube_table = [tube_table; outerjoin(competition_table, chasing_table, ...
            'Keys', 'Mouse_RFID', 'MergeKeys', true)];
    end

    % After all cohorts are processed, assign the concatenated table to the summary_data structure
    summary_data.tube.(range_name) = tube_table;
end

%% Save Updated Summary Data
save(summary_file, 'summary_data');
fprintf('Summary data saved to %s/n', summary_file);
