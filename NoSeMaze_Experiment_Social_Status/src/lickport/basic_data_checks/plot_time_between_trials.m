function [lickport_trials_per_day] = plot_time_between_trials(data_dir,save_dir,cohortsTbl,exclusion_table)
%%

% Predefine empty vectors
timeChar = [];
RewardedTrials = [];
NumLicks = [];
TrialType = [];
ant_lick_matrix = [];

% Loop over cohorts
for i = 1:height(cohortsTbl)
    % Clearing
    clear exclusion_vector d

    % Current cohort
    cohort = cohortsTbl.cohort{i};
    disp(['Processing ' cohort]);
   
    % Define path
    lickport_data_dir = fullfile(data_dir,'interim',cohort,'lickport');

    % Load data
    load(fullfile(lickport_data_dir,['reorganizedData_' cohort '_PhaseCleaned.mat']));
    
    % Exclude animals
    for an=1:length(d)
        
        if any(strcmp(exclusion_table.cohort, cohort) & strcmp(exclusion_table.animalID, d(an).events(1).ID))
            exclusion_vector(an)=true;
        else
            exclusion_vector(an)=false;
        end
    end
    % Exclusion
    d=d(~exclusion_vector);

    for an=1:length(d)
        disp(['Processing animal ' num2str(an)]);

        % Given character array
        timeChar = [timeChar;vertcat(d(an).events.time_StartTrial)];
        RewardedTrials = [RewardedTrials;[d(an).events.reward]'];
        NumLicks = [NumLicks;(cellfun(@numel, {d(an).events.ant_lick}))'];
        TrialType = [TrialType;[d(an).events.trialtype]'];

        % Extract all `ant_lick` values as a cell array
        ant_lick_cells = {d(an).events.ant_lick};

        % Find the maximum number of elements in any row
        max_len = 80;

        % Pad each row with NaN and stack them into a matrix
        clear temp_ant_lick_matrix
        temp_ant_lick_matrix = cellfun(@(x) [x, NaN(1, max_len - numel(x))], ant_lick_cells, 'UniformOutput', false);

        % Convert the padded cell array into a matrix
        ant_lick_matrix = [ant_lick_matrix;vertcat(temp_ant_lick_matrix{:})];
    end
end

% Convert char array to string array (row-wise extraction)
timeStrings = cellstr(timeChar); % Convert char array to cell array of strings

% Convert to duration array
timeDurations = duration(timeStrings, 'InputFormat', 'hh:mm:ss.SSS');

% Compute time differences and convert to numeric seconds
timeDiffs_sec = seconds(diff(timeDurations));

%% histogram of time between trials
figure;
h1 = histogram(timeDiffs_sec(RewardedTrials(1:end-1)==1 & NumLicks(1:end-1)<2),'BinEdges',[0:1:30],'EdgeColor','none'); ax=gca; ax.XLim=[0,30];
hold on;
h2 = histogram(timeDiffs_sec(RewardedTrials(1:end-1)==1 & NumLicks(1:end-1)>=2),'BinEdges',[0:1:30],'EdgeColor','none'); ax=gca; ax.XLim=[0,30];
hold on;
h3 = histogram(timeDiffs_sec(RewardedTrials(1:end-1)==0 & NumLicks(1:end-1)<2),'BinEdges',[0:1:30],'EdgeColor','none'); ax=gca; ax.XLim=[0,30];
hold on;
h4 = histogram(timeDiffs_sec(RewardedTrials(1:end-1)==0 & NumLicks(1:end-1)>=2),'BinEdges',[0:1:30],'EdgeColor','none'); ax=gca; ax.XLim=[0,30];
xlabel('time difference (s)');
ylabel('frequency');
title([num2str(height(cohortsTbl)) ' cohorts']);
xlim([0 30])
set(gcf, 'Units', 'centimeters');
set_fonts();
set(gcf, 'Position', [3 3 9 5]);
legend({'miss','hit','corr. rej.','timeout'},'Interpreter','none');

%%
% save
[annot, srcInfo] = docDataSrc(gcf,fullfile(save_dir),mfilename('fullpath'),logical(1));
exportgraphics(gcf, fullfile(save_dir,'histogram_time_between_trials.pdf'),'ContentType','vector','BackgroundColor','none');

% export source data
% writetable(array2table(lickport_trials_per_day),fullfile(save_dir,'SourceData_histogram_time_between_trials.xlsx'));

close all