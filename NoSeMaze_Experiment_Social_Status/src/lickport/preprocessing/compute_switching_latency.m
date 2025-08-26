function lick_params = compute_switching_latency(lickport_data, lick_params)
%% Switching latency: quantify CS re-learning by calculating latency to switch (in number of trials)
% This function computes the latency to switch behavior after a reversal phase,
% by quantifying how quickly the animal adapts its response to the CS+ (rewarded)
% and CS- (non-rewarded) stimuli. It uses the lick data before and after the reversal phase.
edges = 0.5:3.5;  % Time window for counting licks after odor delivery
reversal_index = lick_params.reversal_index;

% Initialize variables to store latencies and lick counts
cs_plus_switch_latency_at_cs = [];
cs_plus_switch_latency_at_us = [];
cs_minus_switch_latency_at_cs = [];
cs_minus_switch_latency_at_us = [];

% Loop over each phase (after the first phase)
for phase = 2:numel(reversal_index)

    % Initialize arrays for trials and lick counts
    cs_plus_trials_before_current_reversal = [];
    cs_plus_trials_after_current_reversal = [];
    cs_minus_trials_before_current_reversal = [];
    cs_minus_trials_after_current_reversal = [];

    cs_plus_lick_count_before_current_reversal = [];
    cs_plus_lick_count_after_current_reversal = [];
    cs_minus_lick_count_before_current_reversal = [];
    cs_minus_lick_count_after_current_reversal = [];

    % Identify the trials before and after the current reversal phase for CS+ (rewarded) trials
    cs_plus_trials_before_current_reversal = reversal_index(phase-1)-1 + ...
        find([lickport_data(reversal_index(phase-1):reversal_index(phase)-1).reward] == 1, 150, 'last')';
    if phase < numel(reversal_index)
        cs_plus_trials_after_current_reversal = reversal_index(phase)-1 + ...
            find([lickport_data(reversal_index(phase):reversal_index(phase+1)-1).reward] == 1)';
    else
        cs_plus_trials_after_current_reversal = reversal_index(phase)-1 + ...
            find([lickport_data(reversal_index(phase):end).reward] == 1)';
    end
    
    % Identify the trials before and after the current reversal phase for CS+ (rewarded) trials
    cs_minus_trials_before_current_reversal = reversal_index(phase-1) - 1 + ...
        find([lickport_data(reversal_index(phase-1):reversal_index(phase)-1).reward] == 0, 150, 'last')';
    if phase < numel(reversal_index)
        cs_minus_trials_after_current_reversal = reversal_index(phase)-1 + ...
                find([lickport_data(reversal_index(phase):reversal_index(phase+1)-1).reward] == 0)';
    else
        cs_minus_trials_after_current_reversal = reversal_index(phase)-1 + ...
                find([lickport_data(reversal_index(phase):end).reward] == 0)';
    end

    %% Count licks before and after the reversal
    for tr = 1:numel(cs_plus_trials_before_current_reversal)
        cs_plus_lick_count_before_current_reversal = cat(1, cs_plus_lick_count_before_current_reversal, histcounts(lickport_data(cs_plus_trials_before_current_reversal(tr)).licks_aft_od, edges));
    end
    for tr = 1:numel(cs_minus_trials_before_current_reversal)
        cs_minus_lick_count_before_current_reversal = cat(1, cs_minus_lick_count_before_current_reversal, histcounts(lickport_data(cs_minus_trials_before_current_reversal(tr)).licks_aft_od, edges));
    end

    for tr = 1:numel(cs_plus_trials_after_current_reversal)
        cs_plus_lick_count_after_current_reversal = cat(1, cs_plus_lick_count_after_current_reversal, histcounts(lickport_data(cs_plus_trials_after_current_reversal(tr)).licks_aft_od, edges));
    end
    for tr = 1:numel(cs_minus_trials_after_current_reversal)
        cs_minus_lick_count_after_current_reversal = cat(1, cs_minus_lick_count_after_current_reversal, histcounts(lickport_data(cs_minus_trials_after_current_reversal(tr)).licks_aft_od, edges));
    end

    %% More stringent thresholding: calculate latency to switch behavior using moving averages
    clear cs_plus_switch_timepoint_CS cs_plus_switch_timepoint_US cs_minus_switch_timepoint_CS cs_minus_switch_timepoint_US
    cs_plus_switch_timepoint_CS = find(movmean(mean(cs_plus_lick_count_after_current_reversal(:,1:2),2) > 0.7 * mean(cs_plus_lick_count_before_current_reversal(:,1:2), 'all'), [0 9]) > 0.6, 1, 'first');
    if isempty(cs_plus_switch_timepoint_CS)
        cs_plus_switch_timepoint_CS = length(cs_plus_lick_count_after_current_reversal);
    end
    cs_plus_switch_timepoint_US = find(movmean(cs_plus_lick_count_after_current_reversal(:,3) > 0.7 * mean(cs_plus_lick_count_before_current_reversal(:,3), 'all'), [0 9]) > 0.6, 1, 'first');
    if isempty(cs_plus_switch_timepoint_US)
        cs_plus_switch_timepoint_US = length(cs_plus_lick_count_after_current_reversal);
    end
    cs_minus_switch_timepoint_CS = find(movmean(mean(cs_minus_lick_count_after_current_reversal(:,1:2),2) < 0.5 * mean(cs_plus_lick_count_before_current_reversal(:,1:2), 'all'), [0 9]) > 0.6, 1, 'first');
    if isempty(cs_plus_switch_timepoint_CS)
        cs_minus_switch_timepoint_CS = length(cs_minus_lick_count_after_current_reversal);
    end
    cs_minus_switch_timepoint_US = find(movmean(cs_minus_lick_count_after_current_reversal(:,3) < 0.5 * mean(cs_plus_lick_count_before_current_reversal(:,3), 'all'), [0 9]) > 0.6, 1, 'first');
    if isempty(cs_plus_switch_timepoint_CS)
        cs_minus_switch_timepoint_US = length(cs_minus_lick_count_after_current_reversal);
    end
    cs_plus_switch_latency_at_cs = cat(1, cs_plus_switch_latency_at_cs, cs_plus_switch_timepoint_CS);
    cs_plus_switch_latency_at_us = cat(1, cs_plus_switch_latency_at_us, cs_plus_switch_timepoint_US);
    cs_minus_switch_latency_at_cs = cat(1, cs_minus_switch_latency_at_cs, cs_minus_switch_timepoint_CS);
    cs_minus_switch_latency_at_us = cat(1, cs_minus_switch_latency_at_us, cs_minus_switch_timepoint_US);

    %% Parse individual switches to lick_params
    lick_params.(['cs_plus_switch_latency_at_cs_rev', num2str(phase-1)]) = cs_plus_switch_timepoint_CS;
    lick_params.(['cs_plus_switch_latency_at_us_rev', num2str(phase-1)]) = cs_plus_switch_timepoint_US;
    lick_params.(['cs_minus_switch_latency_at_cs_rev', num2str(phase-1)]) = cs_minus_switch_timepoint_CS;
    lick_params.(['cs_minus_switch_latency_at_us_rev', num2str(phase-1)]) = cs_minus_switch_timepoint_US;

    %% Simplified NaN Handling
    lick_params = assignNaNIfEmpty(lick_params, ['cs_plus_switch_latency_at_cs_rev', num2str(phase-1)]);
    lick_params = assignNaNIfEmpty(lick_params, ['cs_plus_switch_latency_at_us_rev', num2str(phase-1)]);
    lick_params = assignNaNIfEmpty(lick_params, ['cs_minus_switch_latency_at_cs_rev', num2str(phase-1)]);
    lick_params = assignNaNIfEmpty(lick_params, ['cs_minus_switch_latency_at_us_rev', num2str(phase-1)]);
end

%% Compute mean latencies for CS+ and CS- switch latencies
lick_params.cs_plus_switch_latency_at_cs_mean = nanmean(cs_plus_switch_latency_at_cs);
lick_params.cs_plus_switch_latency_at_us_mean = nanmean(cs_plus_switch_latency_at_us);
lick_params.cs_minus_switch_latency_at_cs_mean = nanmean(cs_minus_switch_latency_at_cs);
lick_params.cs_minus_switch_latency_at_us_mean = nanmean(cs_minus_switch_latency_at_us);
lick_params.cs_plus_switch_latency_at_cs_median = nanmedian(cs_plus_switch_latency_at_cs);
lick_params.cs_plus_switch_latency_at_us_median = nanmedian(cs_plus_switch_latency_at_us);
lick_params.cs_minus_switch_latency_at_cs_median = nanmedian(cs_minus_switch_latency_at_cs);
lick_params.cs_minus_switch_latency_at_us_median = nanmedian(cs_minus_switch_latency_at_us);

%% Compute CS+ and CS- switch latency shaping
% This function calculates the switch latency shaping for both CS+ and CS- trials
% by computing the ratio of the switch latencies during the first reversal
% to those from the last reversal.
%
% The switch latency shaping is calculated as:
phase = numel(reversal_index);
try
    lick_params.cs_plus_switch_latency_at_cs_shaping = (lick_params.cs_plus_switch_latency_at_cs_rev1 - lick_params.(['cs_plus_switch_latency_at_cs_rev', num2str(phase-1)]))/(lick_params.cs_plus_switch_latency_at_cs_rev1 + lick_params.(['cs_plus_switch_latency_at_cs_rev', num2str(phase-1)]));
    lick_params.cs_plus_switch_latency_at_us_shaping = (lick_params.cs_plus_switch_latency_at_us_rev1 - lick_params.(['cs_plus_switch_latency_at_us_rev', num2str(phase-1)]))/(lick_params.cs_plus_switch_latency_at_us_rev1 + lick_params.(['cs_plus_switch_latency_at_us_rev', num2str(phase-1)]));
    lick_params.cs_minus_switch_latency_at_cs_shaping = (lick_params.cs_minus_switch_latency_at_cs_rev1 - lick_params.(['cs_minus_switch_latency_at_cs_rev', num2str(phase-1)]))/(lick_params.cs_minus_switch_latency_at_cs_rev1 + lick_params.(['cs_minus_switch_latency_at_cs_rev', num2str(phase-1)]));
    lick_params.cs_minus_switch_latency_at_us_shaping = (lick_params.cs_minus_switch_latency_at_us_rev1 - lick_params.(['cs_minus_switch_latency_at_us_rev', num2str(phase-1)]))/(lick_params.cs_minus_switch_latency_at_us_rev1 + lick_params.(['cs_minus_switch_latency_at_us_rev', num2str(phase-1)]));
catch
    % If any required data is missing or calculation fails, set to NaN
    lick_params.cs_plus_switch_latency_at_cs_shaping = NaN;
    lick_params.cs_plus_switch_latency_at_us_shaping = NaN;
    lick_params.cs_minus_switch_latency_at_cs_shaping = NaN;
    lick_params.cs_minus_switch_latency_at_us_shaping = NaN;
end

end

% Simplify NaN assignment by creating a utility function
function lick_params = assignNaNIfEmpty(lick_params, fieldName)
if isempty(lick_params.(fieldName)) || isnan(lick_params.(fieldName))
    lick_params.(fieldName) = NaN;
end
end