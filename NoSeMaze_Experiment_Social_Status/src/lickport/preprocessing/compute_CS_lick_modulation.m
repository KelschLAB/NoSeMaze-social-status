function lick_params = compute_CS_lick_modulation(lickport_data, lick_params)
%% CS+ and CS- lick modulation during the last 150 trials before reversal
% This function calculates the lick modulation for CS+ and CS- based on the lick counts
% during specific windows, such as before and after the reversal.

% Note: Animals might learn sequences in the trial.

% Estimate baseline lick rates for trials before reversal separately for cs
% plus and cs minus trials (baseline lick rates might differ, see
% below)
% Identify reversal indices
reversal_index = [1, find(diff([lickport_data.phase]) == 1) + 1];
lick_params.reversal_index = reversal_index;

% Define trials of interest
% 1. before reversal
cs_plus_trials_before_reversal = [];
cs_minus_trials_before_reversal = [];
% Find trials for CS+ and CS- before reversal
for rev = 2:numel(reversal_index)
    cs_plus_trials_before_reversal = cat(1, cs_plus_trials_before_reversal, reversal_index(rev-1) - 1 + find([lickport_data(reversal_index(rev-1):reversal_index(rev)-1).reward] == 1, 150, 'last')');
    cs_minus_trials_before_reversal = cat(1, cs_minus_trials_before_reversal, reversal_index(rev-1) - 1 + find([lickport_data(reversal_index(rev-1):reversal_index(rev)-1).reward] == 0, 150, 'last')');
end
lick_params.cs_plus_trials_before_reversal = cs_plus_trials_before_reversal;
lick_params.cs_minus_trials_before_reversal = cs_minus_trials_before_reversal;

% baseline_licks in the 150 trials before reversal
baseline_licks = {lickport_data(cs_minus_trials_before_reversal).licks_bef_od};
lick_params.baseline_rate_CSminus_mean = mean(cell2mat(cellfun(@numel, baseline_licks,'UniformOutput',0)).*2);
baseline_licks = cellfun(@(x) x(x>0.05),baseline_licks,'UniformOutput',0);
lick_params.baseline_rate_CSminus_mean_omitfirst = mean(cell2mat(cellfun(@numel, baseline_licks,'UniformOutput',0)).*(1./0.45));% changed by JR: duration is not 500 ms, but 450 ms in here
baseline_licks = {lickport_data(cs_plus_trials_before_reversal).licks_bef_od};
lick_params.baseline_rate_CSplus_mean = mean(cell2mat(cellfun(@numel, baseline_licks,'UniformOutput',0)).*2);
baseline_licks = cellfun(@(x) x(x>0.05),baseline_licks,'UniformOutput',0);
lick_params.baseline_rate_CSplus_mean_omitfirst = mean(cell2mat(cellfun(@numel, baseline_licks,'UniformOutput',0)).*(1./0.45));

%
cs_plus_lick_count = [];
cs_minus_lick_count = [];
cs_plus_lick_count_delta = [];
cs_minus_lick_count_delta = [];
edges = 1:0.05:2.4; %edges(1) = [];

for tr = 1:numel(cs_plus_trials_before_reversal)
    % high-temporal resolution psth over window
    cs_plus_lick_count = cat(1,cs_plus_lick_count, histcounts(lickport_data(cs_plus_trials_before_reversal(tr)).licks_aft_od,edges));
    % trial-by-trial delta to baseline
    cs_plus_lick_count_delta(tr) = nnz(lickport_data(cs_plus_trials_before_reversal(tr)).licks_aft_od>1 & lickport_data(cs_plus_trials_before_reversal(tr)).licks_aft_od<2.4)/1.4-numel(lickport_data(cs_plus_trials_before_reversal(tr)).licks_bef_od>0.05).*(1/0.45);
end
for tr = 1:numel(cs_minus_trials_before_reversal)
    cs_minus_lick_count = cat(1,cs_minus_lick_count, histcounts(lickport_data(cs_minus_trials_before_reversal(tr)).licks_aft_od,edges));
    cs_minus_lick_count_delta(tr) = nnz(lickport_data(cs_minus_trials_before_reversal(tr)).licks_aft_od>1 & lickport_data(cs_minus_trials_before_reversal(tr)).licks_aft_od<2.4)/1.4-numel(lickport_data(cs_minus_trials_before_reversal(tr)).licks_bef_od>0.05).*(1/0.45);
end
lick_params.cs_plus_modulation_averaged_to_base = mean(sum(cs_plus_lick_count,2)./1.4)/lick_params.baseline_rate_CSplus_mean_omitfirst;
lick_params.cs_minus_modulation_averaged_to_base = mean(sum(cs_minus_lick_count,2)./1.4)/lick_params.baseline_rate_CSminus_mean_omitfirst;
lick_params.cs_plus_modulation_trialwise = median(cs_plus_lick_count_delta);
lick_params.cs_minus_modulation_trialwise = median(cs_minus_lick_count_delta);
lick_params.cs_plus_modulation_peak_to_base = max(mean(cs_plus_lick_count,1)*20)/lick_params.baseline_rate_CSplus_mean_omitfirst; % % 20 is correct!
lick_params.cs_minus_modulation_peak_to_base = max(mean(cs_minus_lick_count,1)*20)/lick_params.baseline_rate_CSminus_mean_omitfirst; %
lick_params.cs_plus_modulation_min_to_base = min(mean(cs_plus_lick_count,1)*20)/lick_params.baseline_rate_CSplus_mean_omitfirst; %
lick_params.cs_minus_modulation_min_to_base = min(mean(cs_minus_lick_count,1)*20)/lick_params.baseline_rate_CSminus_mean_omitfirst; %
lick_params.cs_plus_modulation_averaged = mean(sum(cs_plus_lick_count,2));
lick_params.cs_minus_modulation_averaged = mean(sum(cs_minus_lick_count,2));
lick_params.cs_plus_modulation_peak = max(mean(cs_plus_lick_count,1)*20);
lick_params.cs_minus_modulation_peak = max(mean(cs_minus_lick_count,1)*20);
lick_params.cs_plus_modulation_min = min(mean(cs_plus_lick_count,1)*20);
lick_params.cs_minus_modulation_min = min(mean(cs_minus_lick_count,1)*20);
lick_params.cs_plus_modulation_averaged_minus_base = mean(sum(cs_plus_lick_count,2)./1.4)-lick_params.baseline_rate_CSplus_mean_omitfirst;
lick_params.cs_minus_modulation_averaged_minus_base = mean(sum(cs_minus_lick_count,2)./1.4)-lick_params.baseline_rate_CSminus_mean_omitfirst;
lick_params.cs_plus_modulation_peak_minus_base = max(mean(cs_plus_lick_count,1)*20)-lick_params.baseline_rate_CSplus_mean_omitfirst; %
lick_params.cs_minus_modulation_peak_minus_base = max(mean(cs_minus_lick_count,1)*20)-lick_params.baseline_rate_CSminus_mean_omitfirst; %
lick_params.cs_plus_modulation_min_minus_base = min(mean(cs_plus_lick_count,1)*20)-lick_params.baseline_rate_CSplus_mean_omitfirst; %
lick_params.cs_minus_modulation_min_minus_base = min(mean(cs_minus_lick_count,1)*20)-lick_params.baseline_rate_CSminus_mean_omitfirst; %


% CS- lick modulation: Compute the lick modulation during the CS- period, normalized by baseline rate
lick_params.cs_minus_modulation_full_window_averaged = ...
    mean(cs_minus_lick_count, 'all') / lick_params.baseline_rate_CSminus_mean_omitfirst;
end
