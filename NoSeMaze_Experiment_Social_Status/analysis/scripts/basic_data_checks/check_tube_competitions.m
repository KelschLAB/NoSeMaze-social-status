%% Script performs basic data checks/plots, e.g., number of events per day and daycycle

% Pre-Clearing
% clear all

% Specify the main working directory where all the experiment data and configurations are located
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory

% Add necessary paths for the source code. 'genpath' ensures subdirectories under 'src/tube/preprocessing' are included
addpath(genpath(fullfile(main_dir,'src','tube','basic_data_checks')));
addpath(genpath(fullfile(main_dir,'src','helpers')));

% Read cohort information from a configuration file. The file is expected to be a CSV located in the 'config' folder
cohortsTbl = readtable(fullfile(main_dir,'config','cohorts_info.csv'));

% Filter the table to include only cohorts with tube data of interest (use_tube == 1)
% cohortsTbl = cohortsTbl(cohortsTbl.use_tube==1,:);

% Data directory
data_dir = fullfile(main_dir,'data');

% Ouput directory
save_dir = fullfile(main_dir,'results','figures','cross_cohort','basic_data_checks','tube');
if ~isfolder(save_dir)
    mkdir(save_dir)
end

% Plot of tube chasings
[competitions_per_day] = plot_tube_competitions_n_per_day(data_dir,save_dir,cohortsTbl);

% save competition data
writetable([cohortsTbl(:,1),array2table(competitions_per_day)],fullfile(save_dir,'competitions_per_day.csv'));
writetable([cohortsTbl(:,1),array2table([nanmean(competitions_per_day,2),nanmedian(competitions_per_day,2)],'VariableNames',{'mean','median'})],fullfile(save_dir,'competitions_per_day.csv'));
disp(['!!! Check median/mean data in ' save_dir ' and decide which cohorts to use !!!']);

% Plot of tube chasings
plot_tube_competitions_daycycle(data_dir,save_dir,cohortsTbl);


