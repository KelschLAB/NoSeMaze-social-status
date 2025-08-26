function [NaT_sum] = convert_lickport_data(inputPath, outputPath, cohort, displayOn)

% Information:

% "Convert_Lickport_Data" converts NoSeMaze output csv-files in which lickport behavioral data
% is saved (csv-file) into files compatible to data structure established
% by Max and Laurens (mat-file; "events"-variable). Therefore, data saved in

% Input
% inputPath - path belonging to NoSeMaze output file (csv-file)
% outputPath   - directory in which output (mat file) is saved

% Output
% multiple mat-files (one for each animal taking part in current
% experiment)
% "events" variable saved in each mat-file: mat structure including
% all information about trials performed in the NoSeMaze ...



%%
if nargin < 2 % in case there is no input ...
    % define inputPath ...
    inputPath = uigetdir('Please select input directory.');
    % define outputPath ...
    outputPath = uigetdir('Please select output directory.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% GETTING STARTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Read NoSeMaze lickport data and load it into matlab
lickport_files = dir(fullfile(inputPath,'*.csv'));

InputData = table();

% Loop over index files and load them
for i = 1:length(lickport_files)
    current_data = readtable(fullfile(lickport_files(i).folder,lickport_files(i).name));
    InputData = vertcat(InputData,current_data);
end

% Check Input Data for NaT-values
NaT_sum = sum(isnat(InputData.b_timestamp));
if isempty(NaT_sum)
    disp([cohort ': No Not-a-Time-values found in the input data.'])
else
    disp([cohort ': ' num2str(NaT_sum) ' Not-a-Time-values (= trials) found in the input data and deleted.'])
    InputData = InputData(~isnat(InputData.b_timestamp),:);
end

% info: InputData
% column 1 = ID
% column 2 = timestamp
% column 3 = num licks
% ...
% column 5 = reward delivrery + current trial type


%% 2. Get list/number of animals in current experiment (prep animal loop)
% info: ID = 'default' means that an unambiguous assignment was not
% possible in the current trial ...

% create list of all  animals ...
animal_IDs  = unique(InputData.a_animal_id);
% remove 'default' from 'animal_IDs' ...
animal_IDs(strcmp(animal_IDs,'default')) = [];
% get number of animals in current experiment ...
NumAnimalIDs = numel(animal_IDs);


%% 3. Loop over animals
for animal_cur = 1:NumAnimalIDs
    
    % clear variables for the loop run %
    clear ID date time_StartTrial curr_odor_num events

    % get ID of current animal ...
    IDcur = animal_IDs{animal_cur};

    % select all trials performed by current animal ...
    TrialIndex_cur= find(strcmp(InputData.a_animal_id, IDcur));
    TrialNum_cur = numel(TrialIndex_cur);

    % Loop over trials
    for trial_cur = 1:TrialNum_cur
        if displayOn
            disp(['Processing animal ' IDcur ', trial ' num2str(trial_cur)])
        end
        % get row in which current trial data is saved ...
        row_cur = TrialIndex_cur(trial_cur);

        %extract all information of interest from NoSeMaze data and save
        %it in 'events' variable (structure) ...
        % ID
        events(trial_cur).ID = IDcur;

        % date
        events(trial_cur).date = datestr(InputData.b_timestamp(row_cur),'yyyy-mm-dd');

        % time_StartTrial
        events(trial_cur).time_StartTrial= datestr(InputData.b_timestamp(row_cur),'HH:MM:SS.FFF');
        
        % curr_odor_num
        events(trial_cur).curr_odor_num = str2double(InputData.e_rewarded{row_cur}(end));

        % curr_odor_num
        if strcmp(InputData.k_timeout{row_cur},'False')
            events(trial_cur).timeout = 0;
        elseif strcmp(InputData.k_timeout{row_cur},'True')
            events(trial_cur).timeout = 1;
        else
            events(trial_cur).timeout = NaN;
        end
        if strcmp(InputData.h_response_l{row_cur},'False')
            events(trial_cur).response = 0;
        elseif strcmp(InputData.h_response_l{row_cur},'True')
            events(trial_cur).response = 1;
        else
            events(trial_cur).response = NaN;
        end

        % reward_active
        events(trial_cur).reward=str2double(InputData.e_rewarded{row_cur}(1));
        % events(trial_cur).tester(:,1)=str2double(InputData.e_rewarded{row_cur}(1));
        % events(trial_cur).tester(:,2)=str2double(InputData.e_rewarded{row_cur}(5));
        % events(trial_cur).tester(:,3) = str2double(InputData.e_rewarded{row_cur}(end));

        % fv_on
        events(trial_cur).fv_on=0.5;

        % fv_off
        events(trial_cur).fv_off=2.5;

        % rew_delay
        events(trial_cur).rew_delay=2.5;

        %% licks_bef_od
        licks_bef_od=InputData.l_wait_timestamps_l(row_cur);
        % assignment for trials without licking
        if strcmp(licks_bef_od,'not licked')
            events(trial_cur).licks_bef_od=[];
            % assignment for licks:
            % licks are sometimes saves in 'HH:mm:ss.SSS' format (option 1), thus reflecting a time of the day,
            % sometimes as relative time in seconds in relation to the start time, e.g, 0.243|0.398|... (option 2).
        else
            % option 1: 'HH:mm:ss.SSS' format
            if ischar(licks_bef_od{1}) && contains(licks_bef_od{1}, ':')

                % Handle exact timestamp case
                if displayOn
                    disp(['Exact timestamp detected: ' licks_bef_od{1}]);
                end
                % Timestamp as cell
                timestampCell = licks_bef_od;

                % Convert StartTime to a string with the format 'HH:mm:ss.SSS'
                startTimeStr = datestr(InputData.b_timestamp(row_cur), 'HH:MM:SS.FFF');  % 'FFF' for milliseconds

                % Concatenate with a fixed date (ensure a space between date and time)
                StartTime = datetime(['2025-01-01 ' startTimeStr], 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

                % Extract timestamp string from the cell array
                timestampsStr = timestampCell{1};

                % Split the timestamp string into individual timestamps
                timestampsSplit = strsplit(timestampsStr, '|');

                % Concatenate the fixed date to each timestamp
                timestamps = datetime(cellfun(@(x) ['2025-01-01 ' x], timestampsSplit, 'UniformOutput', false), ...
                    'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

                % Check the result of concatenation for debugging
                % disp(timestamps);

                % Calculate the time difference (this will give the result in 'duration' format)
                time_diff = timestamps - StartTime;

                % Convert the 'duration' to seconds and store it as a double
                time_diff_in_seconds = seconds(time_diff);

                % Display the result as seconds, including milliseconds
                % disp(time_diff_in_seconds);

                % Store time differences in events(trial_cur).licks_bef_od
                events(trial_cur).licks_bef_od = time_diff_in_seconds;

                % option 2: time in seconds in relation to the start time, e.g, 0.243|0.398|..
            elseif isnumeric(str2double(licks_bef_od{1}))

                % Handle relative time point case
                if displayOn
                    disp(['Relative time point detected: ' licks_bef_od{1}]);
                end
                % Convert the relative time point to a numeric value
                timestampsSplit = strsplit(licks_bef_od{1}, '|');

                % Store time differences in events(trial_cur).licks_bef_od
                events(trial_cur).licks_bef_od = str2double(timestampsSplit);
            else
                error('Unexpected format in licks_bef_od.');
            end
        end

        %% licks_aft_od
        licks_aft_od=InputData.n_lick_timestamp_l(row_cur);

        % assignment for trials without licking
        if strcmp(licks_aft_od,'not licked')
            events(trial_cur).licks_aft_od=[];

            % assignment for licks:
            % licks are sometimes saves in 'HH:mm:ss.SSS' format (option 1), thus reflecting a time of the day,
            % sometimes as relative time in seconds in relation to the start time, e.g, 0.243|0.398|... (option 2).
        else
            % option 1: 'HH:mm:ss.SSS' format
            if ischar(licks_aft_od{1}) && contains(licks_aft_od{1}, ':')

                % Handle exact timestamp case
                if displayOn
                    disp(['Exact timestamp detected: ' licks_aft_od{1}]);
                end

                % Timestamp as cell
                timestampCell = licks_aft_od;

                % Convert StartTime to a string with the format 'HH:mm:ss.SSS'
                startTimeStr = datestr(InputData.b_timestamp(row_cur), 'HH:MM:SS.FFF');  % 'FFF' for milliseconds

                % Concatenate with a fixed date (ensure a space between date and time)
                StartTime = datetime(['2025-01-01 ' startTimeStr], 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

                % Extract timestamp string from the cell array
                timestampsStr = timestampCell{1};

                % Split the timestamp string into individual timestamps
                timestampsSplit = strsplit(timestampsStr, '|');

                % Concatenate the fixed date to each timestamp
                timestamps = datetime(cellfun(@(x) ['2025-01-01 ' x], timestampsSplit, 'UniformOutput', false), ...
                    'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

                % Check the result of concatenation for debugging
                % disp(timestamps);

                % Calculate the time difference (this will give the result in 'duration' format)
                time_diff = timestamps - StartTime;

                % Convert the 'duration' to seconds and store it as a double
                time_diff_in_seconds = seconds(time_diff);

                % Display the result as seconds, including milliseconds
                % disp(time_diff_in_seconds);

                % Store time differences in events(trial_cur).licks_bef_od
                events(trial_cur).licks_aft_od = time_diff_in_seconds;

                % option 2: time in seconds in relation to the start time, e.g, 0.243|0.398|..
            elseif isnumeric(str2double(licks_aft_od{1}))

                % Handle relative time point case
                if displayOn
                    disp(['Relative time point detected: ' licks_aft_od{1}]);
                end

                % Convert the relative time point to a numeric value
                timestampsSplit = strsplit(licks_aft_od{1}, '|')

                % Store time differences in events(trial_cur).licks_bef_od
                events(trial_cur).licks_aft_od = str2double(timestampsSplit);
            else
                error('Unexpected format in licks_aft_od.');
            end
        end

        % drop_or_not
        ant_licks_count(trial_cur)=length(find(events(trial_cur).licks_aft_od <= events(trial_cur).rew_delay));

        if events(trial_cur).reward==1
            events(trial_cur).drop_or_not='Drop';
        else
            events(trial_cur).drop_or_not='No Drop';
        end

        if ant_licks_count(trial_cur)<=1
            events(trial_cur).drop_or_not='No Drop';
        end

        events(trial_cur).drop_or_not=replace(events(trial_cur).drop_or_not,'No Drop','0');
        events(trial_cur).drop_or_not=replace(events(trial_cur).drop_or_not,'Drop','1');
        events(trial_cur).drop_or_not=str2num(events(trial_cur).drop_or_not);


        % licked_at_the_false_odor

    end % of trial loop

    %%save 'events' in mat-file%%
    savestring_cur = fullfile(outputPath, ['lickport_' cohort '_' IDcur '.mat']);
    save(savestring_cur, 'events');
end % of animal loop




end % function



