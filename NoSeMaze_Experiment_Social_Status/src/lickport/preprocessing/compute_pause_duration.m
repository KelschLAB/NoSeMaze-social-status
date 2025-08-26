function lick_params = compute_pause_duration(lickport_data, lick_params)
    %% Compute_pause_duration
    % This function identifies "giving-up" episodes where the animal shows a
    % drastic reduction in licking behavior during CS+ or US trials after a reversal. 
    % The licking behavior is assessed by calculating lick counts before and after
    % a reversal phase for CS+ (rewarded) and CS- (non-rewarded) trials.
    
    % Get reversal indices (phases where the contingencies reverse)
    reversal_index = lick_params.reversal_index;

    % Define the window for CS+ licking behavior assessment (time window: 1.5s to 3.5s)
    edges = 1.5:3.5;

    % Loop through each reversal phase and compute the licking behavior before and after reversal
    for phase = 2:numel(reversal_index)  % Start from phase 2, so the first-to-last reversal phase
        % Initialize arrays to hold lick counts before and after the reversal
        cs_plus_lick_count_before_current_reversal = [];
        cs_plus_lick_count_after_current_reversal = [];

        % Get CS+ trials before and after the current reversal
        cs_plus_trials_before_current_reversal = reversal_index(phase-1)-1 + ...
            find([lickport_data(reversal_index(phase-1):reversal_index(phase)).reward] == 1, 150, 'last')';
        if phase < numel(reversal_index)
            cs_plus_trials_after_current_reversal = reversal_index(phase)-1 + find([lickport_data(reversal_index(phase):reversal_index(phase+1)).reward] == 1)';
        else
            cs_plus_trials_after_current_reversal = reversal_index(phase)-1 + find([lickport_data(reversal_index(phase):end).reward] == 1)';
        end


        % Count the licks for CS+ trials before the reversal phase
        for tr = 1:numel(cs_plus_trials_before_current_reversal)
            cs_plus_lick_count_before_current_reversal = cat(1, cs_plus_lick_count_before_current_reversal, ...
                histcounts(lickport_data(cs_plus_trials_before_current_reversal(tr)).licks_aft_od, edges));
        end

        % Count the licks for CS+ trials after the reversal phase
        for tr = 1:numel(cs_plus_trials_after_current_reversal)
            cs_plus_lick_count_after_current_reversal = cat(1, cs_plus_lick_count_after_current_reversal, ...
                histcounts(lickport_data(cs_plus_trials_after_current_reversal(tr)).licks_aft_od, edges));
        end

        % Calculate the number of trials where licking during CS+ drops below 25% of pre-reversal frequency
        lick_params.(['pause_duration_at_CS_rev', num2str(phase-1)]) = find(movmean(cs_plus_lick_count_after_current_reversal(:,1) > 0.25*mean(cs_plus_lick_count_before_current_reversal(:,1)), [0 9]) > .6, 1, 'first');
        lick_params.(['pause_duration_at_US_rev', num2str(phase-1)]) = find(movmean(cs_plus_lick_count_after_current_reversal(:,2) > 0.25*mean(cs_plus_lick_count_before_current_reversal(:,2)), [0 9]) > .6, 1, 'first');

        % If no value was found for one of the parameters, set it equal to the other or NaN if neither are available
        if isempty(lick_params.(['pause_duration_at_CS_rev', num2str(phase-1)])) && ~isempty(lick_params.(['pause_duration_at_US_rev', num2str(phase-1)]))
            lick_params.(['pause_duration_at_CS_rev', num2str(phase-1)]) = lick_params.(['pause_duration_at_US_rev', num2str(phase-1)]);
        elseif isempty(lick_params.(['pause_duration_at_CS_rev', num2str(phase-1)])) && isempty(lick_params.(['pause_duration_at_US_rev', num2str(phase-1)]))
            lick_params.(['pause_duration_at_CS_rev', num2str(phase-1)]) = NaN;
        end

        if isempty(lick_params.(['pause_duration_at_US_rev', num2str(phase-1)])) && ~isempty(lick_params.(['pause_duration_at_CS_rev', num2str(phase-1)]))
            lick_params.(['pause_duration_at_US_rev', num2str(phase-1)]) = lick_params.(['pause_duration_at_CS_rev', num2str(phase-1)]);
        elseif isempty(lick_params.(['pause_duration_at_US_rev', num2str(phase-1)])) && ~isempty(lick_params.(['pause_duration_at_CS_rev', num2str(phase-1)]))
            lick_params.(['pause_duration_at_US_rev', num2str(phase-1)]) = NaN;
        end

    end


    % Calculate the ratio of pause durations for CS+ and US at different phases of reversal shaping
    phase = numel(reversal_index);
    % CS
    lick_params.pause_duration_at_CS_shaping_rev_1to2 = 100 * ((lick_params.pause_duration_at_CS_rev1 - lick_params.pause_duration_at_CS_rev2) / ...
        lick_params.pause_duration_at_CS_rev1);
    lick_params.pause_duration_at_CS_shaping_rev_1toLast = 100 * ((lick_params.pause_duration_at_CS_rev1 - lick_params.(['pause_duration_at_CS_rev' num2str(phase-1)])) / ...
        lick_params.pause_duration_at_CS_rev1);
    lick_params.pause_duration_at_CS_shaping_rev_2toLast = 100 * ((lick_params.pause_duration_at_CS_rev2 - lick_params.(['pause_duration_at_CS_rev' num2str(phase-1)])) / ...
        lick_params.pause_duration_at_CS_rev2);
    % US
    lick_params.pause_duration_at_US_shaping_rev_1to2 = 100 * ((lick_params.pause_duration_at_US_rev1 - lick_params.pause_duration_at_US_rev2) / ...
        lick_params.pause_duration_at_US_rev1);
    lick_params.pause_duration_at_US_shaping_rev_1toLast = 100 * ((lick_params.pause_duration_at_US_rev1 - lick_params.(['pause_duration_at_US_rev' num2str(phase-1)])) / ...
        lick_params.pause_duration_at_US_rev1);
    lick_params.pause_duration_at_US_shaping_rev_2toLast = 100 * ((lick_params.pause_duration_at_US_rev2 - lick_params.(['pause_duration_at_US_rev' num2str(phase-1)])) / ...
        lick_params.pause_duration_at_US_rev2);

end
