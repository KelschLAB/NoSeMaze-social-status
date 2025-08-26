function lick_params = compute_giving_up_episodes(lickport_data, lick_params)
%% Identify giving-up episodes
% This function identifies "giving-up" episodes during reversals based on licking behavior.
% It checks if the animal's licking behavior decreases after a reversal at the CS+ (rewarded)
% condition, which could indicate that the animal "gave up" on responding.

% Extract the reversal indices from the lick_params structure
reversal_index = lick_params.reversal_index;

% Define the time window (1.5s to 3.5s after odor onset) to count licks
edges = 1.5:3.5;  % Licking is counted in this time window from 1.5s to 3.5s after odor delivery

% Loop over each phase after the first one, as phase 1 is the initial state
for phase = 2:numel(reversal_index)

    % Initialize arrays for storing trial data
    cs_plus_trials_before_current_reversal = [];
    cs_plus_trials_after_current_reversal = [];
    cs_plus_lick_count_before_current_reversal = [];
    cs_plus_lick_count_after_current_reversal = [];

    % Identify the trials before and after the current reversal phase for CS+ (rewarded) trials
    cs_plus_trials_before_current_reversal = reversal_index(phase-1)-1 + ...
        find([lickport_data(reversal_index(phase-1):reversal_index(phase)-1).reward] == 1, 150, 'last')';
    if phase < numel(reversal_index)
        cs_plus_trials_after_current_reversal = reversal_index(phase)-1 + ...
            find([lickport_data(reversal_index(phase):reversal_index(phase+1)-1).reward] == 1, 70, 'first')';
    else
        cs_plus_trials_after_current_reversal = reversal_index(phase)-1 + ...
            find([lickport_data(reversal_index(phase):end).reward] == 1, 70, 'first')';
    end

    % Remove the first 20 trials after reversal to avoid initial artifacts
    cs_plus_trials_after_current_reversal(1:20) = [];

    % Count the number of licks in the pre-reversal period (before the reversal)
    for tr = 1:numel(cs_plus_trials_before_current_reversal)
        cs_plus_lick_count_before_current_reversal = cat(1, cs_plus_lick_count_before_current_reversal, ...
            histcounts(lickport_data(cs_plus_trials_before_current_reversal(tr)).licks_aft_od, edges));
    end

    % Count the number of licks in the post-reversal period (after the reversal)
    for tr = 1:numel(cs_plus_trials_after_current_reversal)
        cs_plus_lick_count_after_current_reversal = cat(1, cs_plus_lick_count_after_current_reversal, ...
            histcounts(lickport_data(cs_plus_trials_after_current_reversal(tr)).licks_aft_od, edges));
    end

    %% Determine if the animal gave up at the CS (Conditioned Stimulus) and US (Unconditioned Stimulus)
    % Giving up is defined as a significant decrease in licks after reversal compared to before reversal.

    % If the average lick count after reversal is less than 25% of the average lick count before reversal,
    % mark it as a giving-up episode at CS.
    if mean(cs_plus_lick_count_after_current_reversal(:, 1)) < 0.25 * mean(cs_plus_lick_count_before_current_reversal(:, 1))
        lick_params.(['giving_up_at_CS_rev', num2str(phase-1)]) = 1;
    else
        lick_params.(['giving_up_at_CS_rev', num2str(phase-1)]) = 0;
    end

    % Similarly, if the average lick count at US is less than 25% of the average lick count before reversal,
    % mark it as a giving-up episode at US.
    if mean(cs_plus_lick_count_after_current_reversal(:, 2)) < 0.25 * mean(cs_plus_lick_count_before_current_reversal(:, 2))
        lick_params.(['giving_up_at_US_rev', num2str(phase-1)]) = 1;
    else
        lick_params.(['giving_up_at_US_rev', num2str(phase-1)]) = 0;
    end
end
end
