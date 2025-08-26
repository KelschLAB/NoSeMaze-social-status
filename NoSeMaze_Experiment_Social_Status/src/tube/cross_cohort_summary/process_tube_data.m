function metrics_table = process_tube_data(full_hierarchy, day_range, cohortDetails, matrix_type, suffix, plot_dir)
% Processes tube data for a given day range and appends a suffix to metrics.
%
% Parameters:
% full_hierarchy - Hierarchy data structure for the cohort
% day_range - Range of days for the current period
% cohortDetails - Table with basic animal details
% matrix_type - Either 'match_matrix' or 'match_matrix_chasing'
% suffix - String appended to metric names (e.g., '_Competition', '_Chasing')
%
% Returns:
% metrics_table - Table with suffixed metrics aligned to cohortDetails

% Check if day_range is fully covered in full_hierarchy
is_empty = cellfun(@(x) isempty(x.ID), num2cell(full_hierarchy));
empty_days = find(is_empty);

if all(ismember(day_range,empty_days))
    fprintf(['No data available for this day range in full_hierarchy of ' cohortDetails.cohort{1} ', ' suffix(2:end) '. Filling with NaN.']);
    metrics_table = add_empty_metrics(cohortDetails, suffix); % Add empty metrics columns
    return;
end

% Combine match matrices over the day range
current_match_matrix = sum(cat(3, full_hierarchy(day_range).(matrix_type)), 3);

% Get the IDs from the hierarchy data
hierarchy_IDs = full_hierarchy(day_range(end)).ID; % Assuming IDs remain consistent across days in the range

% Compute metrics
metrics = compute_tube_metrics(current_match_matrix,plot_dir);

% Add suffix to metric field names
metrics = add_suffix_to_metrics(metrics, suffix);

% Convert metrics to a table
metrics_table_raw = struct2table(metrics);

% Add RFID column to metrics table
metrics_table_raw.Mouse_RFID = hierarchy_IDs; % Directly from full_hierarchy(day_range(1)).ID


% Handle missing RFIDs (e.g., empty cells)
is_missing = cellfun(@isempty, metrics_table_raw.Mouse_RFID);
metrics_table_raw.Mouse_RFID(is_missing) = {'Unknown'}; % Replace missing RFIDs with 'Unknown'

% Ensure consistent data type for RFIDs
metrics_table_raw.Mouse_RFID = cellfun(@char, metrics_table_raw.Mouse_RFID, 'UniformOutput', false);

% Align metrics with cohortDetails based on Mouse_RFID
metrics_table = outerjoin(cohortDetails, metrics_table_raw, ...
    'Keys', 'Mouse_RFID', 'MergeKeys', true, 'Type', 'Left');

%% Count detections
% Step 1: Find indices where Data_filtered is NOT empty
validIndices = find(~arrayfun(@(i) isempty(full_hierarchy(i).Data_filtered), 1:numel(full_hierarchy)));
% Step 2: Create a new structure with only valid entries
full_hierarchy_temp = full_hierarchy(validIndices);
cellArray = arrayfun(@(i) full_hierarchy_temp(i).Data_filtered.Animal, 1:length(validIndices), 'UniformOutput', false);
all_detections = vertcat(cellArray{:});
all_detections = all_detections(ismember(all_detections,hierarchy_IDs(~cellfun(@isempty, hierarchy_IDs))));

[G_C,G_ID,G_P] = groupcounts(all_detections);
detections_table = table();
detections_table.detection_count = G_C;
detections_table.detection_percentage = G_P;
detections_table.detections_per_day = G_C./length(validIndices);
detections_table.Mouse_RFID = G_ID;

% Check that G_ID and hierarchy_IDs are in the same order
if ~isequal(G_ID, hierarchy_IDs)
    disp('Count IDs and hierarchy IDs are different. Only existing IDs will be merged.');
end

% Align metrics with cohortDetails based on Mouse_RFID
metrics_table = outerjoin(metrics_table, detections_table, ...
    'Keys', 'Mouse_RFID', 'MergeKeys', true, 'Type', 'Left');
end

function metrics_with_suffix = add_suffix_to_metrics(metrics, suffix)
% Adds a suffix to the field names of a metrics struct.
%
% Parameters:
% metrics - Struct containing metric fields (e.g., DSz, Rank)
% suffix - String to append as a suffix to each field name
%
% Returns:
% metrics_with_suffix - Struct with suffixed field names

metrics_with_suffix = struct();
fields = fieldnames(metrics);
for i = 1:numel(fields)
    new_field_name = [fields{i}, suffix]; % Append suffix
    metrics_with_suffix.(new_field_name) = metrics.(fields{i}); % Copy data
end
end

function empty_table = add_empty_metrics(cohortDetails, suffix)
% Adds empty metric columns (filled with NaN) to cohortDetails.
%
% Parameters:
% cohortDetails - Table with basic animal details
% suffix - Suffix to append to metric column names (e.g., '_Competition')
%
% Returns:
% empty_table - Table with added metric columns filled with NaN

% Define metric names
metric_names = {'DS','DSz', 'Rank', 'Fraction_Wins', 'Fraction_Losses', 'N_Wins', 'N_Losses', 'N_Events', 'N_WinsMinusLosses_normalized', ...
    'boxcox_Fraction_Wins', 'boxcox_Fraction_Losses', 'boxcox_N_Wins', 'boxcox_N_Losses', 'boxcox_N_Events',...
    'cuberoot_Fraction_Wins', 'cuberoot_Fraction_Losses', 'cuberoot_N_Wins', 'cuberoot_N_Losses', 'cuberoot_N_Events',...
    'log10_Fraction_Wins', 'log10_Fraction_Losses', 'log10_N_Wins', 'log10_N_Losses', 'log10_N_Events'};
suffixed_names = strcat(metric_names, suffix); % Add suffix to metric names
suffixed_names = [suffixed_names, 'detection_count', 'detection_percentage', 'detections_per_day'];

% Add empty columns to cohortDetails
for i = 1:numel(suffixed_names)
    cohortDetails.(suffixed_names{i}) = NaN(height(cohortDetails), 1); % Fill with NaN
end

empty_table = cohortDetails; % Return the updated table
end
