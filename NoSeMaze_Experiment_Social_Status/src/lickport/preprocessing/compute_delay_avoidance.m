function lick_params = compute_delay_avoidance(lickport_data, lick_params)
    %% Delay Avoidance Learner: Determine if the animal has learned to avoid the delay
    % This function calculates whether the animal has learned to avoid licking
    % during CS- trials (non-rewarded trials), which is interpreted as the
    % animal successfully learning the delay avoidance task.
    %
    % It checks if the animal has made fewer than or equal to 1 lick, on average, 
    % across all CS- trials during the defined period after odor delivery.
    % 
    % Inputs:
    % - lickport_data: A structure with lick data for each trial.
    % - lick_params: A structure that stores the calculated learning parameters.
    %
    % Outputs:
    % - lick_params: Updated structure with the delay avoidance learning status.

    % Define the edges for counting licks within a 0.5s to 2.5s window
    edges = 0.5:2.5;  % Time bins for counting licks after odor delivery

    %% 1. Calculation of delay avoidance in CS- trials over the whole period
    % Extract CS- trials before the reversal from lick_params
    cs_minus_trials_before_reversal = lick_params.cs_minus_trials_before_reversal;

    % Initialize an array to store lick counts for CS- trials before reversal
    cs_minus_lick_count_before = [];

    % Count licks for each CS- trial before reversal
    for tr = 1:numel(cs_minus_trials_before_reversal)
        % Count the licks within the defined time window (edges)
        cs_minus_lick_count_before = cat(1, cs_minus_lick_count_before, ...
            histcounts(lickport_data(cs_minus_trials_before_reversal(tr)).licks_aft_od, edges));
    end

    % Check if the average lick count is low (less than or equal to 1) across the CS- trials
    % If the animal refrains from licking, it is considered as having learned the task
    dal = mean(cs_minus_lick_count_before, 'all') <= 1;  % Delay avoidance learner check

    % Update lick_params with the result of the delay avoidance learning check
    if dal
        lick_params.delay_avoidance_learner = 1;  % Success: Animal has learned delay avoidance
    else
        lick_params.delay_avoidance_learner = 0;  % Failure: Animal has not learned delay avoidance
    end

    
    %% 2. Calculation of delay avoidance in CS- trials before the first and the last reversal
    % Count licks for each CS- trial before reversal in the first and the
    % last phase

    % Initialize arrays for trials and lick counts
    cs_minus_trials_before_first_reversal = [];
    cs_minus_lick_count_before_first_reversal = [];
    cs_minus_trials_before_last_reversal = [];
    cs_minus_lick_count_before_last_reversal = [];

    % Identify the trials before and after the current reversal phase for CS-(non-rewarded) trials
    % First reversal (phase = 2)
    reversal_index = lick_params.reversal_index;
    phase = 2;
    cs_minus_trials_before_first_reversal = reversal_index(phase-1) - 1 + ...
        find([lickport_data(reversal_index(phase-1):reversal_index(phase)-1).reward] == 0, 150, 'last')';
    % Last reversal
    phase = numel(reversal_index);
    cs_minus_trials_before_last_reversal = reversal_index(phase-1) - 1 + ...
        find([lickport_data(reversal_index(phase-1):reversal_index(phase)-1).reward] == 0, 150, 'last')';

    % Count licks for each CS- trial before first reversal
    for tr = 1:numel(cs_minus_trials_before_first_reversal)
        % Count the licks within the defined time window (edges)
        cs_minus_lick_count_before_first_reversal = cat(1, cs_minus_lick_count_before_first_reversal, ...
            histcounts(lickport_data(cs_minus_trials_before_first_reversal(tr)).licks_aft_od, edges));
    end

       % Count licks for each CS- trial before last reversal
    for tr = 1:numel(cs_minus_trials_before_last_reversal)
        % Count the licks within the defined time window (edges)
        cs_minus_lick_count_before_last_reversal = cat(1, cs_minus_lick_count_before_last_reversal, ...
            histcounts(lickport_data(cs_minus_trials_before_last_reversal(tr)).licks_aft_od, edges));
    end

    % Check if the average lick count is low (less than or equal to 1) across the CS- trials
    % If the animal refrains from licking, it is considered as having learned the task
    dal_first = mean(cs_minus_lick_count_before_first_reversal, 'all') <= 1;  % Delay avoidance learner check
    dal_last = mean(cs_minus_lick_count_before_last_reversal, 'all') <= 1;  % Delay avoidance learner check

    % Update lick_params with the result of the delay avoidance learning check
    if dal_first
        lick_params.delay_avoidance_learner_firstReversal = 1;  % Success: Animal has learned delay avoidance
    else
        lick_params.delay_avoidance_learner_firstReversal = 0;  % Failure: Animal has not learned delay avoidance
    end
    if dal_last
        lick_params.delay_avoidance_learner_lastReversal = 1;  % Success: Animal has learned delay avoidance
    else
        lick_params.delay_avoidance_learner_lastReversal = 0;  % Failure: Animal has not learned delay avoidance
    end

    % Compute the ratio of lick behavior before the first reversal and
    % before the last reversal
    lick_params.delay_avoidance_shaping = (mean(cs_minus_lick_count_before_first_reversal, 'all') - mean(cs_minus_lick_count_before_last_reversal, 'all')) / ...
        (mean(cs_minus_lick_count_before_first_reversal, 'all') + mean(cs_minus_lick_count_before_last_reversal, 'all'));

    % If the result is infinite (e.g., division by zero), set it to NaN
    if isinf(lick_params.delay_avoidance_shaping)
        lick_params.delay_avoidance_shaping = NaN;
    end
end
