%% main_preprocessing_lickport_data.m

% Pre-Clearing
clear all

% Specify the main working directory where all the experiment data and configurations are located
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory

% Add necessary paths for the source code. 'genpath' ensures subdirectories under 'src/lickport/preprocessing' are included
addpath(genpath(fullfile(main_dir,'src','lickport','preprocessing')));
addpath(genpath(fullfile(main_dir,'src','helpers')));

% Read cohort information from a configuration file. The file is expected to be a CSV located in the 'config' folder
cohortsTbl = readtable(fullfile(main_dir,'config','cohorts_info.csv'));

% Animal exclusion list
% Sometimes animals have to be excluded, e.g., because they only have very
% few trials or were only recognized during a part of the time.
exclusion_table = table({'cohort09','cohort11','cohort11','cohort12','cohort13','cohort14','cohort15','cohort17','cohort18','cohort19','cohort20'}',{'0007CDEAA7','0007CB0B74','0007F2B3FF','0007CB0AD4','0007F2E455','0007F2F283','0007CB0B74','0007F2F283','0007CB0AD4','0007CB0AD4','0007CB0DBC'}','VariableNames',{'cohort','animalID'});
% cohort09, 0007CDEAA7: only phase 1 and 2
% cohort11, 0007CB0B74: less than 120 trials
% cohort11, 0007F2B3FF: less than 1500 trials, only one phase
% cohort12, 0007CB0AD4: only 505 trials
% cohort13, 0007F2E455: less than 2000 trials
% cohort14, 0007F2F283: less than 200 trials
% cohort15, 0007CB0B74: less than 50 trials
% cohort17, 0007F2F283: less than 200 trials
% cohort18, 0007CB0AD4: very few trials in the first weeks
% cohort19, 0007CB0AD4: only 200 trials
% cohort20, 0007CB0DBC: only trials for phase 1 and phase 2

% Predefine lick parameters for all cohorts
lick_params = [];

% Loop through each cohort listed in the cohorts table
for i = 1:height(cohortsTbl)
    % Set time
    tic

    % Extract the current cohort's name from the table
    cohort = cohortsTbl.cohort{i};

    % Display current cohort 
    disp(['Processing ' cohort])

    % Define paths for raw and interim data for this cohort
    rawPath = fullfile(main_dir,'data','raw',cohort);
    interimPath = fullfile(main_dir,'data','interim',cohort);
    processedPath = fullfile(main_dir,'data','processed',cohort);

    %% Preprocess lickport data
    % Begin the preprocessing pipeline for lickport data

    %% 1. Convert csv file to mat-data for each animal

    % convert_lickport_data converts raw NoSeMaze CSV files containing
    % lickport behavioral data from the NoSeMaze into lickport_<cohort>_<animal_ID>.mat files.
    % The converted data adheres to the format established by Max and Laurens,
    % with a structured events variable that contains trial data for each animal.
    if 1==0

        % Define paths for raw lickport data and the location to save processed LOG files
        rawLickportPath = fullfile(rawPath,'lickport');
        interimLickportPath = fullfile(interimPath,'lickport');
        if ~isfolder(interimLickportPath)
            mkdir(interimLickportPath);
        end

        % Display option (e.g., shows animal ID, trial, data format that is processed)
        displayOn = false; % true

        % Convert raw CSV files to LOG format
        [NaT_sum] = convert_lickport_data(rawLickportPath,interimLickportPath,cohort,displayOn);

        % Missing values
        MissingValues(i) = NaT_sum;
    end

    %% 2. Combine lickport data, add additional information (e.g.,
    % anticipatory licks), save in the interimLickportPath a
    % "reorganizedData_cohortXX.mat"-file
    if 1==0
        % Path where interim lickport data is stored
        interimLickportPath = fullfile(interimPath,'lickport');
        combine_lickport_data_single_animals(interimLickportPath,cohort);
    end

    %% 3. Clean lickport data
    if 1==0
        % Path where interim lickport data is stored
        interimLickportPath = fullfile(interimPath, 'lickport');  % Constructing the full path to the interim data directory

        % Get the file information of the reorganized data specific to the cohort
        d_struc_file = dir(fullfile(interimLickportPath, ['reorganizedData_' cohort '.mat']));

        % Load the reorganized data file for the specified cohort
        load(fullfile(d_struc_file.folder, d_struc_file.name));

        % Call a custom function to clean phases and exclude certain animals from the data
        % The function `clean_phases_and_exclude_animals` is modifying `d` (the dataset),
        % cleaning phases and excluding animals based on the provided cohort and exclusion table
        [d] = clean_phases_and_exclude_animals(d, cohort, exclusion_table);

        % Save the cleaned data into a new file with '_PhaseCleaned' appended to the filename  
        save(fullfile(d_struc_file.folder, ['reorganizedData_' cohort '_PhaseCleaned.mat']), 'd');
    end

    %% 4. Plot PSTH 
    if 1==0
        % Path where interim lickport data is stored
        interimLickportPath = fullfile(interimPath,'lickport');
        % Path for plots
        plotPath = fullfile(main_dir,'results','figures',cohort,'basic_data_plots','lickport','plots_PSTH');
        if ~isfolder(plotPath)
            mkdir(plotPath);
        end
        % Load d structure file
        d_struc_file = dir(fullfile(interimLickportPath,['reorganizedData_' cohort '_PhaseCleaned.mat']));
        load(fullfile(d_struc_file.folder,d_struc_file.name));
        % Plots
        plot_PSTH_single_animals(d, plotPath, cohort);
    end
                
    %% 5. Compute lick features
    % Compute lick features 
    if 1==1
        % Path where interim lickport data is stored
        interimLickportPath = fullfile(interimPath,'lickport');
        % Path for plots
        plotPath = fullfile(interimLickportPath,'plots_baseline_over_task');
        if ~isfolder(plotPath)
            mkdir(plotPath);
        end
        % Load d structure file
        d_struc_file = dir(fullfile(interimLickportPath,['reorganizedData_' cohort '_PhaseCleaned.mat']));
        load(fullfile(d_struc_file.folder,d_struc_file.name));
        % 
        lick_params_current_cohort = compute_learning_parameters(d,cohort,plotPath);
        
        % Find missing fields in lick_params
        if exist('lick_params', 'var') && ~isempty(lick_params)
            missing_fields = setdiff(fieldnames(lick_params), fieldnames(lick_params_current_cohort));

            % Add missing fields to lick_params with default NaN values
            for i = 1:length(missing_fields)
                [lick_params_current_cohort(:).(missing_fields{i})] = deal(NaN);  % Assign NaN to all existing elements
            end

            missing_fields = setdiff(fieldnames(lick_params_current_cohort), fieldnames(lick_params));

            % Add missing fields to lick_params with default NaN values
            for i = 1:length(missing_fields)
                [lick_params(:).(missing_fields{i})] = deal(NaN);  % Assign NaN to all existing elements
            end
        end

        % Update lick parameters
        lick_params = [lick_params,lick_params_current_cohort];

        % Save lick parameters
        if ~isfolder(fullfile(processedPath,'lickport'))
            mkdir(fullfile(processedPath,'lickport'));
        end
        save(fullfile(processedPath,'lickport',['lick_params_' cohort '.mat']),'lick_params_current_cohort')
    end

    toc
end

% Save lick parameters
outputDir = fullfile(main_dir,'data','processed','cross_cohort_files','lickport');
if ~isfolder(outputDir)
    mkdir(outputDir);
end
save(fullfile(outputDir,'lick_params_cross_cohort.mat'),'lick_params')
