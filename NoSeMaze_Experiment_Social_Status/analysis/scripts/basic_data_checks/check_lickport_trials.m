%% Script performs basic data checks/plots, e.g., number of lick events per day and daycycle

% Pre-Clearing
clear all

% Specify the main working directory where all the experiment data and configurations are located
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory

% Add necessary paths for the source code. 'genpath' ensures subdirectories under 'src/tube/preprocessing' are included
addpath(genpath(fullfile(main_dir,'src','lickport','basic_data_checks')));
addpath(genpath(fullfile(main_dir,'src','helpers')));

% Read cohort information from a configuration file. The file is expected to be a CSV located in the 'config' folder
cohortsTbl = readtable(fullfile(main_dir,'config','cohorts_info.csv'));

% Filter the table to include only cohorts with tube data of interest (use_tube == 1)
cohortsTbl = cohortsTbl(cohortsTbl.use_lickport==1,:);

% Data directory
data_dir = fullfile(main_dir,'data');

% Ouput directory
save_dir = fullfile(main_dir,'results','figures','cross_cohort','basic_data_checks','lickport');
if ~isfolder(save_dir)
    mkdir(save_dir)
end

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

% Plot of trials per day
if 1==0
    [lickport_trials_per_day] = plot_lickport_trials_n_per_day(data_dir,save_dir,cohortsTbl,exclusion_table);
end

% Plot of trials over daycycle
if 1==1
    plot_lickport_trials_daycycle(data_dir,save_dir,cohortsTbl,exclusion_table);
end

% Plot pauses between trials
if 1==0
    plot_time_between_trials(data_dir,save_dir,cohortsTbl,exclusion_table)
end
