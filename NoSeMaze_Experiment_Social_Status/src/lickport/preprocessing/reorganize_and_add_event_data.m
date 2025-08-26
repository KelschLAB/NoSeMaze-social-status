function [d] = reorganize_and_add_event_data(inputData, outputDir, cohort)
    % GNG_REORD_FUN - Reorganizes the event data, adjusts lick data, and adds phase and trial type information.
    % 
    % Inputs:
    %   inputData - The input data (structure array) containing events data.
    %   outputDir - Directory path to save the processed data.
    %   cohort   
    % 
    % Outputs:
    %   d - The reorganized data with added trial types and phases.
    
    % Step 1: Field renaming and reordering
    d = inputData;
    
    % Step 2: Adjust lick data (both pre- and post-reward)
    for i = 1:numel(d)
        for j = 1:numel(d(i).events)
            if isempty(d(i).events(j).licks_aft_od)
                d(i).events(j).ant_lick = [];
                d(i).events(j).post_lick = [];
            else
                % Separate licks before and after the reward delay
                d(i).events(j).ant_lick = d(i).events(j).licks_aft_od(d(i).events(j).licks_aft_od <= d(i).events(j).rew_delay);
                d(i).events(j).post_lick = d(i).events(j).licks_aft_od(d(i).events(j).licks_aft_od > d(i).events(j).rew_delay);
            end
        end

        for j = 1:numel(d(i).events)
            if isempty(d(i).events(j).licks_aft_od)
                d(i).events(j).ant_lick_count = 0;
                d(i).events(j).post_lick_count = 0;
            else
                d(i).events(j).ant_lick_count = sum(d(i).events(j).licks_aft_od <= d(i).events(j).rew_delay);
                d(i).events(j).post_lick_count = sum(d(i).events(j).licks_aft_od > d(i).events(j).rew_delay);
            end
        end
    end
    
    % Step 3: Assign trial types based on lick counts and rewards
    for i = 1:numel(d)
        for j = 1:numel(d(i).events)
            if d(i).events(j).reward == 1
                if d(i).events(j).ant_lick_count >= 2
                    d(i).events(j).trialtype = 1;  % Reward, and at least 2 anticipatory licks
                else
                    d(i).events(j).trialtype = 2;  % Reward, but less than 2 anticipatory licks
                end
            else
                if d(i).events(j).ant_lick_count >= 2
                    d(i).events(j).trialtype = 4;  % No reward, but at least 2 anticipatory licks
                else
                    d(i).events(j).trialtype = 3;  % No reward, and less than 2 anticipatory licks
                end
            end
        end
    end
    
    % Step 4: Phase specification based on odor and reward conditions
    for i = 1:numel(d)
        if d(i).events(1).reward == 1
            odorReward = d(i).events(1).curr_odor_num;
            odorNoReward = setdiff([1, 2], odorReward);  % Other odor
        else
            odorNoReward = d(i).events(1).curr_odor_num;
            odorReward = setdiff([1, 2], odorNoReward);  % Other odor
        end
        
        % Assign phase based on odor and reward conditions
        aux_chap = zeros(1, numel(d(i).events));
        for j = 1:numel(d(i).events)
            if d(i).events(j).curr_odor_num == odorReward && d(i).events(j).reward == 1
                aux_chap(j) = 1;  % Reward odor
            elseif d(i).events(j).curr_odor_num == odorNoReward && d(i).events(j).reward == 0
                aux_chap(j) = 1;  % No reward odor
            else
                aux_chap(j) = 2;  % Other cases (non-reward)
            end
        end
        
        % Define phases based on changes in odor-reward conditions
        phase = ones(size(aux_chap));  % Default phase value
        changeIndices = find(diff(aux_chap) ~= 0) + 1;  % Indices where phase changes
        for k = 1:length(changeIndices)
            phase(changeIndices(k):end) = k + 1;  % Increment phase value
        end
        
        % Assign phase to events
        for j = 1:numel(d(i).events)
            d(i).events(j).phase = phase(j);
            d(i).events(j).ph_rev = aux_chap(j);
        end
    end
    
    % Step 5: Save the reorganized data to a specified file path
    filename = ['reorganizedData_' cohort '.mat'];
    if ~isfolder(outputDir)
        mkdir(outputDir);  % Create directory if it doesn't exist
    end
    save(fullfile(outputDir, filename), 'd');
end
