function lick_params = compute_time_to_criterion(lickport_data, lick_params)
%% Time to criterion: quantify CS re-learning by calculating the time until animals reach an 80% criterion (in number of trials)
% This function computes the time animals need until they reach a criterion of 80% correct hits/rejections. 
edges = 0.5:2.5;  % Time window for counting licks after odor delivery
reversal_index = lick_params.reversal_index;

cs_plus_time_to_criterion = [];
cs_minus_time_to_criterion = [];

% Loop over each phase
for phase = 1:numel(reversal_index)

    % Initialize arrays for trials and lick counts
    cs_plus_trials_current_phase = [];
    cs_minus_trials_current_phase = [];

    cs_plus_lick_count_current_phase = [];
    cs_minus_lick_count_current_phase = [];

    % Identify the trials after the current phase for CS+ (rewarded) trials
    if phase < numel(reversal_index)
        cs_plus_trials_current_phase = reversal_index(phase)-1 + ...
            find([lickport_data(reversal_index(phase):reversal_index(phase+1)).reward] == 1)';
    else
        cs_plus_trials_current_phase = reversal_index(phase)-1 + ...
            find([lickport_data(reversal_index(phase):end).reward] == 1)';
    end
    
    % Identify the trials before and after the current reversal phase for CS- (non-rewarded) trials
    if phase < numel(reversal_index)
        cs_minus_trials_current_phase = reversal_index(phase)-1 + ...
            find([lickport_data(reversal_index(phase):reversal_index(phase+1)).reward] == 0)';
    else
        cs_minus_trials_current_phase = reversal_index(phase)-1 + ...
            find([lickport_data(reversal_index(phase):end).reward] == 0)';
    end

    %% Count licks in the current phase
    for tr = 1:numel(cs_plus_trials_current_phase)
        cs_plus_lick_count_current_phase = cat(1, cs_plus_lick_count_current_phase, histcounts(lickport_data(cs_plus_trials_current_phase(tr)).licks_aft_od, edges));
    end
    for tr = 1:numel(cs_minus_trials_current_phase)
        cs_minus_lick_count_current_phase = cat(1, cs_minus_lick_count_current_phase, histcounts(lickport_data(cs_minus_trials_current_phase(tr)).licks_aft_od, edges));
    end

   %% Calculate time to criterion using moving averages
   clear current_cs_plus_time_to_criterion_timepoint current_cs_minus_time_to_criterion_timepoint
   current_cs_plus_time_to_criterion_timepoint = find(movmean(sum(cs_plus_lick_count_current_phase(:,1:2),2)>=2,[0,9])>0.8,1,'first');
   current_cs_minus_time_to_criterion_timepoint = find(movmean(sum(cs_minus_lick_count_current_phase(:,1:2),2)<2,[0,9])>0.8,1,'first');

   if isempty(current_cs_plus_time_to_criterion_timepoint)
       cs_plus_time_to_criterion = length(cs_plus_trials_current_phase);
       lick_params.(['cs_plus_time_to_criterion_phase', num2str(phase)]) = length(cs_plus_trials_current_phase);
   else
       cs_plus_time_to_criterion = cat(1, cs_plus_time_to_criterion, current_cs_plus_time_to_criterion_timepoint);
       lick_params.(['cs_plus_time_to_criterion_phase', num2str(phase)]) = current_cs_plus_time_to_criterion_timepoint;
   end
   if isempty(current_cs_minus_time_to_criterion_timepoint)
       cs_minus_time_to_criterion = length(cs_minus_trials_current_phase);
       lick_params.(['cs_minus_time_to_criterion_phase', num2str(phase)]) = length(cs_minus_trials_current_phase);
   else
       cs_minus_time_to_criterion = cat(1, cs_minus_time_to_criterion, current_cs_minus_time_to_criterion_timepoint);
       lick_params.(['cs_minus_time_to_criterion_phase', num2str(phase)]) = current_cs_minus_time_to_criterion_timepoint;
   end
end

%% Compute mean latencies for CS+ and CS- switch latencies
lick_params.cs_plus_time_to_criterion_mean = nanmean(cs_plus_time_to_criterion);
lick_params.cs_minus_time_to_criterion_mean = nanmean(cs_minus_time_to_criterion);
lick_params.cs_plus_time_to_criterion_median = nanmedian(cs_plus_time_to_criterion);
lick_params.cs_minus_time_to_criterion_median = nanmedian(cs_minus_time_to_criterion);
end

