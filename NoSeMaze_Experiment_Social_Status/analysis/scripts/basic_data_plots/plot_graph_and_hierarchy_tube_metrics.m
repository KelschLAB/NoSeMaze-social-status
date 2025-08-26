%% plot_graph_and_hierarchy_tube_metrics.m
% Author: Jonathan Reinwald
% Date: January 2025
% Purpose: Generate and save plots of David's Scores and hierarchy graphs
% for tube competition and chasing metrics across cohorts and day ranges.

%% Pre-Clearing
clear; clc; close all;

%% Set Directories
% Define main and subdirectories for data and results
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory
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
% cohortsTbl = cohortsTbl(cohortsTbl.use_tube == 1, :);

%% Define Parameters
% Specify data types, metrics, and day ranges to analyze
data_types = {'tube'}; % Only analyze tube metrics
metrics.tube = {'Fraction_Wins_Chasing', 'Fraction_Losses_Chasing'};%DS_Competition','DS_Chasing', 'ELO_rand_Competition', 'ELO_nonrand_Competition', ...
                 % 'ELO_rand_Chasing', 'ELO_nonrand_Chasing', 
day_ranges = {'D1_21'};%'D1_14', 
day_ranges_idx = {[1, 21]}; % [1, 14], Corresponding indices for day ranges

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
        for metric_idx = 1:numel(metrics.tube)
            current_metric = metrics.tube{metric_idx};
            
            % Iterate through cohorts
            for cohort_idx = 1:height(cohortsTbl)
                cohort = cohortsTbl.cohort{cohort_idx};
                
                % Create output directory for the cohort
                output_dir = fullfile(main_dir, 'results', 'figures', cohort, 'basic_data_plots', data_type);
                if ~isfolder(output_dir)
                    mkdir(output_dir);
                end
                
                % Extract data and IDs for the current cohort
                current_IDs = summary_data.(data_type).(range_name).Mouse_RFID(...
                    strcmp(summary_data.(data_type).(range_name).cohort, cohort));
                current_data = summary_data.(data_type).(range_name).(current_metric)(...
                    strcmp(summary_data.(data_type).(range_name).cohort, cohort));

                % Skip if all data values are NaN
                if all(isnan(current_data))
                    continue;
                end

                % Plot and save hierarchy score and rank data
                [f, sorted_IDs, sorted_data, sorted_ranks] = plot_hierarchy_scores_in_cohort(...
                    current_data(~isnan(current_data)), ...
                    current_IDs(~isnan(current_data)), cohort, current_metric);

                % Save sorted data and ranks as a CSV
                DataTable = table(sorted_IDs, sorted_data, sorted_ranks, ...
                                  'VariableNames', {'IDs', current_metric, 'rank'});
                writetable(DataTable, fullfile(output_dir, ...
                    ['SourceData_' current_metric '_' cohort '_' range_name '.csv']));

                % Save the hierarchy score plot as a PDF% Save the heatmap
                [annot, srcInfo] = docDataSrc(f,fullfile(output_dir),mfilename('fullpath'),logical(1));
                exportgraphics(f, fullfile(output_dir, ...
                    [current_metric '_' cohort '_' range_name '.pdf']), ...
                    'ContentType', 'vector', 'BackgroundColor', 'none');
                
                % Load and process hierarchy data
                cohort_data_dir = fullfile(processed_dir, cohort, data_type, 'full_hierarchy_files');
                full_hierarchy_file = fullfile(cohort_data_dir, ['full_hierarchy_' cohort '.mat']);
                load(full_hierarchy_file, 'full_hierarchy');

                % Select the appropriate match matrix
                if contains(current_metric, 'Competition')
                    current_match_matrix = sum(cat(3, ...
                        full_hierarchy([day_ranges_idx{day_idx}(1):day_ranges_idx{day_idx}(2)]).match_matrix), 3);
                elseif contains(current_metric, 'Chasing')
                    current_match_matrix = sum(cat(3, ...
                        full_hierarchy([day_ranges_idx{day_idx}(1):day_ranges_idx{day_idx}(2)]).match_matrix_chasing), 3);
                end

                % Select IDs for non-empty rows/columns in the match matrix
                ID_selection = ~cellfun(@isempty, full_hierarchy(day_ranges_idx{day_idx}(2)).ID);
                IDs = full_hierarchy(day_ranges_idx{day_idx}(2)).ID(ID_selection);

                % Plot and save hierarchy graph
                f = plot_hierarchy_graph_in_cohort(...
                    current_match_matrix(ID_selection, ID_selection), ...
                    current_data(~isnan(current_data)), IDs, cohort, current_metric);

                % Save match matrix with IDs as a CSV
                csvData = [[{'ID/IDs'}, IDs']; IDs, ...
                           num2cell(current_match_matrix(ID_selection, ID_selection))];
                writecell(csvData, fullfile(output_dir, ...
                    ['SourceData_GraphPlot_' current_metric(strfind(current_metric,'_')+1:end) '_' cohort '_' range_name '.csv']));

                % Save hierarchy graph plot as a PDF
                [annot, srcInfo] = docDataSrc(f,fullfile(output_dir),mfilename('fullpath'),logical(1));
                exportgraphics(f, fullfile(output_dir, ...
                    ['GraphPlot_' current_metric(strfind(current_metric,'_')+1:end) '_' cohort '_' range_name '.pdf']), ...
                    'ContentType', 'vector', 'BackgroundColor', 'none');
                
                % Close all open figures to save memory
                close all;
            end
        end
    end
end





