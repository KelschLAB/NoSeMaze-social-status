%% process_manual_tube_competition.m
% Author: Jonathan Reinwald
% Date: January 2025
% Purpose:
% - Process manual tube competition (MTT) data from raw Excel files.
% - Compute metrics (David's Score, ranks, wins/losses) and update the 
%   summary data for specified day ranges.

%% Pre-Clearing
clear; clc; close all;

%% Set Directories
% Define main and subdirectories
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory
raw_data_dir = fullfile(main_dir, 'data', 'raw');
processed_dir = fullfile(main_dir, 'data', 'processed');
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');

% Add required paths for helper functions and preprocessing scripts
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));
addpath(genpath(fullfile(main_dir, 'src', 'tube', 'preprocessing')));

% Define day ranges of interest
day_ranges = {'D1_14', 'D1_21'};

%% Load Summary Data
% Ensure summary data file exists, and load it
if ~isfile(summary_file)
    error('Summary file not found. Please run `main_build_cross_cohort_summary` first.');
end
load(summary_file, 'summary_data');

%% Locate Manual Tube Test (MTT) Files
% Recursively search for MTT Excel files in raw data directory
MTT_filelist = dir(fullfile(raw_data_dir, '**', 'MTT_*.xlsx'));

%% Process Each MTT File
for m_idx = 1:numel(MTT_filelist)
    
    % Import MTT data and calculate David's Score
    [MTT_match_matrix, MTT_IDs] = import_MTT(fullfile(MTT_filelist(m_idx).folder, MTT_filelist(m_idx).name));
    MTT_DS = compute_DS_from_match_matrix(MTT_match_matrix).DS;
    [sortedDS, sortingIdx] = sort(MTT_DS, 'descend'); % Sort DS for rank computation
    [~, MTT_rank] = sort(sortingIdx); % Determine ranks based on sorted indices
    MTT_DSz = zscore(MTT_DS);
    MTT_tube_competition_wins = sum(MTT_match_matrix, 2); % Row sums (wins)
    MTT_tube_competition_losses = sum(MTT_match_matrix, 1); % Column sums (losses)

    % Extract cohort name from the file path using regular expression
    MTT_cohort = regexp(MTT_filelist(m_idx).folder, 'cohort\d+', 'match', 'once');
    if isempty(MTT_cohort)
        warning('Cohort name not found for file: %s', MTT_filelist(m_idx).name);
        continue;
    end
    
    % Display the processed cohort for confirmation
    disp(['Processing cohort: ', MTT_cohort]);

    %% Initialize Columns in Summary Data
    % Define the columns for MTT metrics
    colMTT = {'DS_MTT', 'DSz_MTT', 'Rank_MTT', 'N_Wins_MTT', 'N_Losses_MTT', 'Fraction_Wins_MTT', 'Fraction_Losses_MTT'};
    
    for day_idx = 1:numel(day_ranges)
        % Ensure the necessary columns exist in the summary data for the current day range
        for col_idx = 1:numel(colMTT)
            if ~ismember(colMTT{col_idx}, summary_data.tube.(day_ranges{day_idx}).Properties.VariableNames)
                summary_data.tube.(day_ranges{day_idx}).(colMTT{col_idx}) = NaN(height(summary_data.tube.(day_ranges{day_idx})), 1);
            end
        end

        %% Update Summary Data with MTT Metrics
        for j = 1:height(MTT_IDs)
            % Match rows in summary data based on Mouse ID and cohort
            mouseIdx = strcmp(summary_data.tube.(day_ranges{day_idx}).Mouse_RFID, MTT_IDs{j}) & ...
                       strcmp(summary_data.tube.(day_ranges{day_idx}).cohort, MTT_cohort);
            
            % If a match is found, update the fields
            if any(mouseIdx)
                summary_data.tube.(day_ranges{day_idx}).DS_MTT(mouseIdx) = MTT_DS(j);
                summary_data.tube.(day_ranges{day_idx}).DSz_MTT(mouseIdx) = MTT_DSz(j);
                summary_data.tube.(day_ranges{day_idx}).Rank_MTT(mouseIdx) = MTT_rank(j);
                summary_data.tube.(day_ranges{day_idx}).N_Wins_MTT(mouseIdx) = MTT_tube_competition_wins(j);
                summary_data.tube.(day_ranges{day_idx}).N_Losses_MTT(mouseIdx) = MTT_tube_competition_losses(j);
                summary_data.tube.(day_ranges{day_idx}).Fraction_Wins_MTT(mouseIdx) = ...
                    MTT_tube_competition_wins(j) / sum(MTT_tube_competition_wins);
                summary_data.tube.(day_ranges{day_idx}).Fraction_Losses_MTT(mouseIdx) = ...
                    MTT_tube_competition_losses(j) / sum(MTT_tube_competition_losses);
            else
                % Issue a warning if no match is found
                warning('No matching row for ID %s and cohort %s in summary_data.', ...
                        MTT_IDs{j}, MTT_cohort);
            end
        end
    end
end

save(summary_file, 'summary_data');
