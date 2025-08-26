function d = combine_lickport_data_single_animals(inputDataFolder, cohort)
    % combine_lickport_data_single_animals - Combines .mat files from single animal and adds phase information.
    % 
    % Inputs:
    %   inputDataFolder - Path to the folder containing the .mat files.
    %   cohort - Cohort identifier or name to be displayed in the error messages.
    % 
    % Outputs:
    %   d - A structure containing the reorganized data. It is saved in the
    %   inputDataFolder as "reorganizedData_cohort01.mat" (d.mat-file)
    %
    % This function assumes that each .mat file already contains a structure called 'events'.
    % The function reads and combines these files into one unified structure.
    
    % Get a list of all .mat files in the folder that start with 'lickport_'
    matFiles = dir(fullfile(inputDataFolder, 'lickport_*.mat'));

    % Initialize an empty structure array to store the events data from each file
    combinedData = [];  % We initialize this to an empty array, and append data in the loop
    
    % Loop through each .mat file in the folder and load the 'events' structure
    for i = 1:numel(matFiles)
        % Load the data from the current file
        fileData = load(fullfile(matFiles(i).folder,matFiles(i).name));
        
        % Append the 'events' structure from the current file to 'combinedData'
        % This approach can be slow for large datasets, consider preallocating if performance is an issue
        combinedData = [combinedData; fileData];  % Append events data
    end

    %  Data Cleaning: minimum lick-timestamp difference = 50ms
    [combinedData] = clean_lick_data(combinedData);

    % Now that we have combined data, we call another function to reorganize the structure further
    % This is adding more details to the data structure
    % Note that saving of the data in the inputDataFolder is happening
    % here.
    d = reorganize_and_add_event_data(combinedData, inputDataFolder, cohort);

    % The following section performs a quality check and counts errors
    % This checks if certain conditions are met and logs errors based on those conditions
    for i = 1:numel(d)
        % Initialize an error count for each entry in the structure
        errorCount(i) = 0;
        
        % Loop through all 'events' for the current entry
        for j = 1:numel(d(i).events)
            % Condition for detecting certain types of errors (based on specific field values)
            if d(i).events(j).ant_lick_count < 2 && d(i).events(j).reward == 1 && d(i).events(j).drop_or_not == 1
                errorCount(i) = errorCount(i) + 1;
                
                % If an error is detected, record the event index
                if errorCount(i) ~= 0
                    errorIndices{i}(errorCount(i), 1) = j;
                end
            end
        end

        % Display the number of errors found for the current animal (from 'cohort') 
        % This is useful for tracking the data quality for each animal
        disp([num2str(errorCount(i)) ' errors found for ' cohort ' animal # ' d(i).events(1).ID])
    end
end
