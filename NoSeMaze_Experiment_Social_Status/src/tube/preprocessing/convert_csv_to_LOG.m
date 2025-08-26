function convert_csv_to_LOG(inputPath,outputPath,ID_replacement_file,start_date,end_date)
% convert csv-file to LOG-data
% Jonathan Reinwald, 07.01.2024

% Find csv-files
csvFiles = dir(fullfile(inputPath, '*.csv')); % Find all .csv files

% Filter out files containing 'replacement' in their name
csvFiles = csvFiles(~contains({csvFiles.name}, 'replacement', 'IgnoreCase', true));

% Display the filtered list (optional)
disp({csvFiles.name});

% Set a counter
counter = 1;
% Loop over csv-files
for ix = 1:size(csvFiles,1)
    clear tempData
    % Load data
    tempData = readtable(fullfile(csvFiles(ix).folder,csvFiles(ix).name));
    
    for iix = 1:size(tempData,1)
        % Bring into LOG structure
        trial_datetime = char(tempData.Date_time(iix));
        Day{counter} = trial_datetime(1:10);
        Time_hour(counter) = hours(duration(trial_datetime(12:end)));
        Time_sec(counter) = seconds(duration(trial_datetime(12:end)));
        Detector(counter) = tempData.UnitNumber(iix);
        Animal(counter) = tempData.TransponderCode(iix);
        % Update counter
        counter = counter+1;
    end
end

Data_full = table(Day',Time_hour',Time_sec',Detector',Animal','VariableNames',{'Day','Time_hour','Time_sec','Detector','Animal'});

% ID replacement and switch, if necessary:
if ~isempty(ID_replacement_file)
    
    % Load replacement file (all columns are read as strings)
    opts = detectImportOptions(fullfile(ID_replacement_file.folder,ID_replacement_file.name), 'TextType', 'string');
    opts = setvartype(opts, 'string'); % Force all columns to be read as strings
    replacement_table = readtable(fullfile(ID_replacement_file.folder,ID_replacement_file.name), opts);
       
    % Replacement of IDs
    for ix = 1:height(replacement_table)
        % Only perform replacement, if row is not missing
        if ~ismissing(replacement_table.IDs_to_replace(ix))
            
            clear logical_replacement_array
            clear number_replacement_array
            
            logical_replacement_array = ismember(Data_full.Animal,replacement_table.IDs_to_replace(ix));
            number_replacement_array=find(logical_replacement_array);
            for kx = 1:length(number_replacement_array)
                Data_full.Animal(number_replacement_array(kx)) = cellstr(replacement_table.IDs_for_replacement(ix));
            end
        end
    end
    
    % Deletion of IDs
    for ix = 1:height(replacement_table)
        
        if ~ismissing(replacement_table.IDs_to_delete(ix))
            
            clear logical_deletion_array
            
            logical_deletion_array = ismember(Data_full.Animal,replacement_table.IDs_to_delete(ix));
            Data_full(logical_deletion_array,:)=[];
        end
    end
end


% Split into single-day files and save in interim data
unique_days = unique(Data_full.Day); % Days with actual data
all_days = datestr((datetime(start_date):datetime(end_date))', 'yyyy-mm-dd'); % All expected days
missing_days = setdiff(all_days, unique_days); % Find missing days

% Create LOG-files directory
logPath = fullfile(outputPath,'LOG-files');
if ~isfolder(logPath)
    mkdir(logPath);
end

% Process days with data
for dd = 1:numel(unique_days)
    clear Data
    Data = Data_full(contains(Data_full.Day, unique_days{dd}), :);
    save(fullfile(logPath, ['LOG_', unique_days{dd}]), 'Data');
end

% Process missing days
if ~isempty(missing_days)
    fprintf('Missing data for %d days:\n', numel(missing_days));
    disp(missing_days'); % Display the missing days
    
    for dd = 1:numel(missing_days)
        clear Data
        Data = table(string.empty(0, 1), [], [], [], string.empty(0, 1), ...
            'VariableNames', {'Day', 'Time_hour', 'Time_sec', 'Detector', 'Animal'});
        save(fullfile(logPath, ['LOG_', missing_days{dd}]), 'Data');
        fprintf('Created placeholder LOG file for missing day: %s\n', missing_days{dd});
    end
else
    fprintf('No missing data. All expected days are accounted for.\n');
end