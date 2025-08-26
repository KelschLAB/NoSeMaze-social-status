%% main_preprocessing_tube_data.m

% Pre-Clearing
clear all

% Specify the main working directory where all the experiment data and configurations are located
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory

% Add necessary paths for the source code. 'genpath' ensures subdirectories under 'src/tube/preprocessing' are included
addpath(genpath(fullfile(main_dir,'src','tube','preprocessing')));
addpath(genpath(fullfile(main_dir,'src','helpers')));

% Read cohort information from a configuration file. The file is expected to be a CSV located in the 'config' folder
cohortsTbl = readtable(fullfile(main_dir,'config','cohorts_info.csv'));

% Loop through each cohort listed in the cohorts table
for i = 9:height(cohortsTbl)
    % Extract the current cohort's name from the table
    cohort = cohortsTbl.cohort{i};
    disp(['Processing: ' cohort])
    % Define paths for raw and interim data for this cohort
    rawPath = fullfile(main_dir,'data','raw',cohort);
    interimPath = fullfile(main_dir,'data','interim',cohort);
    processedPath = fullfile(main_dir,'data','processed',cohort);
    
    %% Preprocess tube data
    % Begin the preprocessing pipeline for tube data

    %%
    % 1. Convert csv file to LOG data
    % Info: This step includes an option for deletion or merging of
    % specific animals based on their IDs (e.g., when IDs switched during
    % the experiment). The file is in csv format and should be located in 
    % each animal's raw data folder (e.g., /home/jonathan.reinwald/NoSeMaze_Experiment/data/raw/cohort17/tube/ID_replacement_cohort17.csv)
    % containing three columns: IDs_to_delete, IDs_to_replace, IDs_for_replacement 
    if 1==0
        % Define paths for raw tube data and the location to save processed LOG files
        rawTubePath = fullfile(rawPath,'tube');
        interimTubePath = fullfile(interimPath,'tube');
        ID_replacement_file = dir(fullfile(rawTubePath,'ID_replacement*.csv'))
        
        % Convert raw CSV files to LOG format
        convert_csv_to_LOG(rawTubePath,interimTubePath,ID_replacement_file,cohortsTbl.start_date(i),cohortsTbl.end_date(i));
    end
    
    %%
    % 2. Extract tube competitions and plot tube competition events
    if 1==1
        % Path where interim LOG data is stored
        interimTubePath = fullfile(interimPath,'tube');
        
        % List all LOG files in the directory (assuming filenames match the pattern 'LOG*.mat')
        LOG_file_list = dir(fullfile(interimTubePath,'LOG-files','LOG*.mat'));
        
        % Loop through each LOG file to process tube competition data
        for file_num = 1:numel(LOG_file_list)
            disp(['Processing: ' LOG_file_list(file_num).name])
            % Full path to the current LOG file
            input_file = fullfile(LOG_file_list(file_num).folder, LOG_file_list(file_num).name);

            % Load data
            load(input_file,'Data');
            
            % Only perform the extraction of tube competitions and plotting
            % for non-empty data-files 
            % The check if ~isempty(Data) ensures that tube competitions and 
            % plots are only generated for non-empty LOG_XXX.mat files. 
            % This is important because placeholder files were created for missing days.
            
            if ~isempty(Data)
                
                % Extract competition data from the LOG file
                double_detection_threshold = 0.5;
                hierarchy_data = extract_tube_competitions_from_LOG_clean(input_file,double_detection_threshold);
                
                % Generate a unique save directory for plots and processed data, based on the current file name
                [~,cur_name,~] = fileparts(input_file);
                hierarchy_data.save_dir = fullfile(interimTubePath,'plots_tube_competitions',cur_name);
                if ~isfolder(hierarchy_data.save_dir)
                    mkdir(hierarchy_data.save_dir);
                end
                
                % Plot tube competition events based on the extracted hierarchy data
                % plot_tube_competitions(hierarchy_data);
                
                % Save the processed data (hierarchy and DS) to a MAT file in the save directory
                % save(fullfile(hierarchy_data.save_dir,'Data_NEW.mat'),'hierarchy_data');
                
                % Close all open figures to prevent memory issues during the next iteration
                close all
            elseif isempty(Data)
                fprintf('Skipping %s: No data available.\n', input_file);
                continue;
            end
        end
    end
    
    %%
    % 3. Check events and label/note incorrect events
    % !!! Before step 4., label exclude_events.mat/include_events.mat/inverse_events.mat manually
    % based on the plots and save them in the plot-folder !!!
    
    % 4. Combine days of curated hierarchy and exclude/include/inverse the
    % labeled data points
    if 1==0
        % Full path to plots_tube_competitions-path, in which the Data.mat
        % with the DS and the hierarchy_info are saved
        inputPath = fullfile(interimPath,'tube','plots_tube_competitions');
        
        % Output directory for full_hierarchy.mat
        interimOutputPath = fullfile(interimPath,'tube','full_hierarchy_files');
        if ~isfolder(interimOutputPath)
                mkdir(interimOutputPath); 
        end

        % Combine days for hierarchy and thereby exclude/inverse/include
        % the labeled events
        combine_days_for_hierarchy(inputPath,interimOutputPath,cohort,cohortsTbl.start_date(i),cohortsTbl.end_date(i));
    end
    
    % 5. Extract chasing from (the filtered data from the) full_hierarchy-files
    if 1==0
        % Input directory for full_hierarchy.mat
        clear full_hierarchy
        interimInputPath = fullfile(interimPath,'tube','full_hierarchy_files');
        full_hierarchy_file = dir(fullfile(interimInputPath,'full_hierarchy*.mat'));
        load(fullfile(full_hierarchy_file.folder,full_hierarchy_file.name));

        % Specify additional input
        ops.plotter = 1;
        ops.save_dir = fullfile(interimPath,'tube','plots_tube_chasings');
        if ~isfolder(ops.save_dir)
            mkdir(ops.save_dir);
        end
        ops.threshold_lag_at_detector = 1.5; % in seconds
        ops.threshold_lag_through_tube = 2; % in seconds
        ops.start_date = cohortsTbl.start_date(i);
        ops.end_date = cohortsTbl.end_date(i);
        
        % Calculate and add tube chasing events to full_hierarchy
        [full_hierarchy]=extract_tube_chasing_from_full_hierarchy(full_hierarchy,ops);
        
        % Save
        processedOutputPath = fullfile(processedPath,'tube','full_hierarchy_files');
        if ~isfolder(processedOutputPath)
                mkdir(processedOutputPath); 
        end
        save(fullfile(interimInputPath,['full_hierarchy_' cohort '.mat']),'full_hierarchy');
        save(fullfile(processedOutputPath,['full_hierarchy_' cohort '.mat']),'full_hierarchy');
        % 
        close all
    end
    
    % 6. Create sequence lists for R
    if 1==0
        % Input directory for full_hierarchy.mat
        clear full_hierarchy
        interimInputPath = fullfile(interimPath,'tube','full_hierarchy_files');
        full_hierarchy_file = dir(fullfile(interimInputPath,['full_hierarchy_' cohort '.mat']));
        day_range = 'full';
        
        % Output directory
        processedOutputPath = fullfile(processedPath,'tube','sequence_files','tube_competitions');
        if ~isfolder(processedOutputPath)
            mkdir(processedOutputPath);
        end
        
        % Create sequence list for tube competitions
        create_sequence_list_tube_competitions(fullfile(full_hierarchy_file.folder,full_hierarchy_file.name),processedOutputPath,cohort,day_range);
        
        % Output directory
        processedOutputPath = fullfile(processedPath,'tube','sequence_files','tube_chasings');
        if ~isfolder(processedOutputPath)
            mkdir(processedOutputPath);
        end
        
        % Create sequence list for tube competitions
        create_sequence_list_tube_chasings(fullfile(full_hierarchy_file.folder,full_hierarchy_file.name),processedOutputPath,cohort,day_range);
        
    end
end

