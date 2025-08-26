% Pre-Clearing
clear; clc;

% Define k
% k = 'K100';
k = 'Koptimal';

% Define directories
mainDir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory
csvDataDir = fullfile(mainDir, '/data/interim/cross_cohort_files/tube');
summaryDataDir = fullfile(mainDir, '/data/processed/cross_cohort_files/');

% Load summary_data.mat
summaryFile = fullfile(summaryDataDir, 'summary_data.mat');
if exist(summaryFile, 'file')
    load(summaryFile, 'summary_data');
else
    error('summary_data.mat not found in the specified directory.');
end

% List all CohortCharacteristics CSV files
csvFiles = dir(fullfile(csvDataDir, ['CohortCharacteristics_*_' k '_combined.csv']));

% Process each CSV file
for i = 1:length(csvFiles)
    % Extract filename
    fileName = csvFiles(i).name;
    fullPath = fullfile(csvDataDir, fileName);

    % Extract modality and day range from the filename
    modality = extractBetween(fileName, 'CohortCharacteristics_Tube', '_');
    modality = modality{1}(1:end-1); % Convert to lowercase (e.g., 'chasings', 'competitions')
    if contains(fileName, '_D')
        dayRange = extractBetween(fileName, '_D', '_K');
        dayRange = ['D' dayRange{1}];
    elseif contains(fileName, '_Last')
        if contains(fileName, 'Chasings')
            dayRange = extractBetween(fileName, 'Chasings_', '_K');
        elseif contains(fileName, 'Competitions')
            dayRange = extractBetween(fileName, 'Competitions_', '_K');
        end
        dayRange = dayRange{1};
    else
        % Handle the case where neither '_D' nor '_Last' is found in fileName
        disp('Neither ''_D'' nor ''_Last'' found in fileName.');
    end

    % Read CSV file
    data = readtable(fullPath);

    % Check if cohort and ID columns match in summary_data
    if ~isfield(summary_data.tube, dayRange)
        % Initialize structure for this day range if it doesn't exist
        summary_data.tube.(dayRange) = table(summary_data.tube.(dayRange).Mouse_RFID, summary_data.tube.(dayRange).cohort, ...
            'VariableNames', {'Mouse_RFID', 'cohort'});
    elseif ~isfield(summary_data.tube, dayRange)
        % Initialize structure for this day range if it doesn't exist
        summary_data.tube.(dayRange) = table(summary_data.tube.(dayRange).Mouse_RFID, summary_data.tube.(dayRange).cohort, ...
            'VariableNames', {'Mouse_RFID', 'cohort'});
    end

    % Initialize columns with NaN values
    colRandELO = ['ELO_rand_', modality, '_', k];
    colNonRandELO = ['ELO_nonrand_', modality, '_', k];
    if ~ismember(colRandELO, summary_data.tube.(dayRange).Properties.VariableNames)
        summary_data.tube.(dayRange).(colRandELO) = NaN(height(summary_data.tube.(dayRange)), 1);
    end
    if ~ismember(colNonRandELO, summary_data.tube.(dayRange).Properties.VariableNames)
        summary_data.tube.(dayRange).(colNonRandELO) = NaN(height(summary_data.tube.(dayRange)), 1);
    end

    % Match rows in summary_data based on ID (Mouse_RFID) and cohort
    for j = 1:height(data)
        % Find matching row in summary_data
        mouseIdx = strcmp(summary_data.tube.(dayRange).Mouse_RFID, data.ID{j}) & ...
            strcmp(summary_data.tube.(dayRange).cohort, data.cohort{j});

        % If a match is found, update the fields
        if any(mouseIdx)
            summary_data.tube.(dayRange).(colRandELO)(mouseIdx) = data.sorted_randELOscore(j);
            summary_data.tube.(dayRange).(colNonRandELO)(mouseIdx) = data.sorted_non_randELOscore(j);
        else
            warning('No matching row for ID %s and cohort %s in summary_data.', ...
                data.ID{j}, data.cohort{j});
        end

        % Internal control
        if any(mouseIdx)
            if (abs(summary_data.tube.(dayRange).(['DS_', modality])(mouseIdx)) - abs(data.sorted_DS(j)))<0.001
                disp('DS from R and matlab data are matching.')
            else
                error('DS from R and Matlab are not matching.')
            end

        end
    end
end

% Save the updated summary_data.mat
save(summaryFile, 'summary_data');
disp('summary_data.mat has been updated successfully.');
