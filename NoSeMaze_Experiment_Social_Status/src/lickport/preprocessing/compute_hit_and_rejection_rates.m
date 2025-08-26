function lick_params = compute_hit_and_rejection_rates(lickport_data,lick_params)

% Extract the reversal indices from the lick_params structure
reversal_index = lick_params.reversal_index;

% Define the time window (1.5s to 3.5s after odor onset) to count licks
edges = 0.5:2.5;  % Licking is counted in this time window from 1.5s to 3.5s after odor delivery

% Initialize arrays for trials and lick counts over all reversals
cs_plus_lick_count_before_reversal = [];
cs_minus_lick_count_before_reversal = [];

% Loop over each phase after the first one, as phase 1 is the initial state
for phase = 2:numel(reversal_index)

    % Initialize arrays for trials and lick counts
    cs_plus_trials_before_current_reversal = [];
    cs_minus_trials_before_current_reversal = [];

    cs_plus_lick_count_before_current_reversal = [];
    cs_minus_lick_count_before_current_reversal = [];

    % Identify the trials before and after the current reversal phase for CS+ (rewarded) trials
    cs_plus_trials_before_current_reversal = reversal_index(phase-1)-1 + ...
        find([lickport_data(reversal_index(phase-1):reversal_index(phase)-1).reward] == 1, 150, 'last')';
        
    % Identify the trials before and after the current reversal phase for CS- (non rewarded) trials
    cs_minus_trials_before_current_reversal = reversal_index(phase-1) - 1 + ...
        find([lickport_data(reversal_index(phase-1):reversal_index(phase)-1).reward] == 0, 150, 'last')';

    %% Count licks before and after the reversal
    for tr = 1:numel(cs_plus_trials_before_current_reversal)
        cs_plus_lick_count_before_current_reversal = cat(1, cs_plus_lick_count_before_current_reversal, histcounts(lickport_data(cs_plus_trials_before_current_reversal(tr)).licks_aft_od, edges));
    end
    for tr = 1:numel(cs_minus_trials_before_current_reversal)
        cs_minus_lick_count_before_current_reversal = cat(1, cs_minus_lick_count_before_current_reversal, histcounts(lickport_data(cs_minus_trials_before_current_reversal(tr)).licks_aft_od, edges));
    end
    
    %% Parse individual rates to lick_params
    lick_params.(['correct_hit_rate_rev', num2str(phase-1)]) = sum(sum(cs_plus_lick_count_before_current_reversal,2)>=2)./size(cs_plus_lick_count_before_current_reversal,1);
    lick_params.(['correct_rejection_rate_rev', num2str(phase-1)]) = sum(sum(cs_minus_lick_count_before_current_reversal,2)<2)./size(cs_minus_lick_count_before_current_reversal,1);

    %% Update overall data
    cs_plus_lick_count_before_reversal = cat(1, cs_plus_lick_count_before_reversal, cs_plus_lick_count_before_current_reversal);
    cs_minus_lick_count_before_reversal = cat(1, cs_minus_lick_count_before_reversal, cs_minus_lick_count_before_current_reversal);
end

%% all trials before switch (stable phase)
lick_params.correct_hit_rate = sum(sum(cs_plus_lick_count_before_reversal,2)>=2)./size(cs_plus_lick_count_before_reversal,1);
lick_params.correct_rejection_rate = sum(sum(cs_minus_lick_count_before_reversal,2)<2)./size(cs_minus_lick_count_before_reversal,1);

%% all trials
cs_plus_lick_count_all = [];
cs_minus_lick_count_all = [];
cs_plus_trials = find([lickport_data.reward]==1);
cs_minus_trials = find([lickport_data.reward]==0);
%% Count licks before and after the reversal
for tr = 1:numel(cs_plus_trials)
    cs_plus_lick_count_all = cat(1, cs_plus_lick_count_all, histcounts(lickport_data(cs_plus_trials(tr)).licks_aft_od, edges));
end
for tr = 1:numel(cs_minus_trials)
    cs_minus_lick_count_all = cat(1, cs_minus_lick_count_all, histcounts(lickport_data(cs_minus_trials(tr)).licks_aft_od, edges));
end
lick_params.correct_hit_rate_allTrials = sum(sum(cs_plus_lick_count_all,2)>=2)./size(cs_plus_lick_count_all,1);
lick_params.correct_rejection_rate_allTrials = sum(sum(cs_minus_lick_count_all,2)<2)./size(cs_minus_lick_count_all,1);