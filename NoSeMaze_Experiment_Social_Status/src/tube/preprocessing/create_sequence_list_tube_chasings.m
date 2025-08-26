function create_sequence_list_tube_chasings(InputFile,OutputPath,cohort,day_range)
% Jonathan Reinwald, 01/2025

% Info:
% This script creates a sequence list of all chasing events in the order of
% occurence

% Load the full_hierarchy.mat file
load(InputFile);

% Define days of interest
if day_range=='full'
    day_range = 1:length(full_hierarchy);
else
    day_range = [day_range(1):day_range(2)];
end

% Display the name of the loaded file
fprintf('Loaded file: %s\n', InputFile);

% Initialize empty vectors to store winner_exit_times, winner_IDs and
% loser_IDs and days
winner_exit_times = [];
winner_IDs = [];
loser_IDs = [];
date = [];

% Loop through the days of interest
for current_day = day_range
    
    % Check if data for this day is available
    if isempty(full_hierarchy(current_day).Data)
        disp([cohort ': No data for day ' num2str(current_day)])
        continue
    else
        % Find non-empty indices
        clear indices
        indices = find(full_hierarchy(current_day).match_matrix_chasing>0);

        % Initialize empty vectors to store winner_exit_times, winner_IDs and
        % loser_IDs and days (This is needed for the sorting)
        winner_exit_times_currentday = [];
        winner_IDs_currentday = [];
        loser_IDs_currentday = [];

        % Loop through each index
        for idx = 1:length(indices)

            % Check if the cell contains a structure
            if ~isempty(full_hierarchy(current_day).match_info_chasing{indices(idx)}) && isstruct(full_hierarchy(current_day).match_info_chasing{indices(idx)})
                % Extract winner_exit_time and add to the vector
                winner_exit_time = [full_hierarchy(current_day).match_info_chasing{indices(idx)}.winner_exit_time];
                winner_exit_times_currentday = [winner_exit_times_currentday; winner_exit_time'];

                winner_ID = {full_hierarchy(current_day).match_info_chasing{indices(idx)}.winner_ID}';
                winner_IDs_currentday = vertcat(winner_IDs_currentday, winner_ID);

                loser_ID = {full_hierarchy(current_day).match_info_chasing{indices(idx)}.loser_ID}';
                loser_IDs_currentday = vertcat(loser_IDs_currentday, loser_ID);
            else
                % Handle cases where the cell is empty or does not contain the expected structure
                warning('Cell (%d, %d) is empty or does not contain a structure.', r, c);
            end
        end

        % Sorting of the data by the time of the interaction (on the
        % respective day)
        clear winner_exit_times_sorted sorting_idx DateCell
        [winner_exit_times_sorted,sorting_idx] = sort(winner_exit_times_currentday,'ascend');
        winner_exit_times = [winner_exit_times; winner_exit_times_sorted];
        winner_IDs = [winner_IDs; winner_IDs_currentday(sorting_idx)];
        loser_IDs = [loser_IDs; loser_IDs_currentday(sorting_idx)];
        current_date = deblank(full_hierarchy(current_day).Data.Day{1});

        % Given timestamp in seconds
        timestamp_seconds = winner_exit_times_sorted;
        % Reference date
        % Define the input date strings
        if length(current_date) == 8
            % Define the input formats
            input_format = 'yyyyMMdd';
        elseif length(current_date) == 10
            input_format = 'yyyy-MM-dd';
        end

        % Convert to datetime objects
        reference_date = datetime(current_date, 'InputFormat', input_format);
        % Convert seconds to days
        timestamp_days = timestamp_seconds / (24 * 60 * 60);
        % Add the days to the reference date
        resulting_date = repmat(reference_date,size(timestamp_seconds,1),1) + days(timestamp_days);
        % Format the resulting date as YYYY-MM-DD
        formatted_date = datestr(resulting_date, 'yyyy-mm-dd HH:MM:SS');

        date = [date;formatted_date];
    end
end

% Creation of a table for csv-files
clear myTable
myTable = table(winner_IDs,loser_IDs,winner_exit_times,date,'VariableNames',{'winners','losers','daytime','day'});
writetable(myTable,fullfile(OutputPath,['SequenceList_',cohort,'_days' num2str(day_range(1)) 'to' num2str(day_range(end)) '_TubeChasings.csv']));


