function lick_params = compute_baseline_rate(lickport_data,lick_params)
    %% Compute baseline lick rate (Hz)
    % This function computes the baseline lick rate during the baseline period
    % (before odor onset), including mean, standard deviation, median, and quartile values.
    
    % Case 1: Perform statistical analysis on the lick data during the baseline
    % period from 0 to 0.5s after trial-start
    lick_params.baseline_rate_isnormal = 1 - adtest(cell2mat(cellfun(@numel, {lickport_data(:).licks_bef_od}, 'UniformOutput', 0)) * 2);
    lick_params.baseline_rate_mean = mean(cell2mat(cellfun(@numel, {lickport_data(:).licks_bef_od}, 'UniformOutput', 0)) * 2);
    lick_params.baseline_rate_sd = std(cell2mat(cellfun(@numel, {lickport_data(:).licks_bef_od}, 'UniformOutput', 0)) * 2);
    lick_params.baseline_rate_median = median(cell2mat(cellfun(@numel, {lickport_data(:).licks_bef_od}, 'UniformOutput', 0)) * 2);
    lick_params.baseline_rate_quartile = quantile(cell2mat(cellfun(@numel, {lickport_data(:).licks_bef_od}, 'UniformOutput', 0)) * 2, [0.25, 0.75]);

    % Case 2: Average lick-rate (Hz) during baseline from 0.05 to 0.5s after trial-start
    % --> 0.05 s as start point to avoid initial lick peak
    baseline_licks_omitfirst = cellfun(@(x) x(x>0.05),{lickport_data(:).licks_bef_od},'UniformOutput',0);
    lick_params.baseline_rate_mean_omitfirst = mean(cell2mat(cellfun(@numel, baseline_licks_omitfirst,'UniformOutput',0)).*(1./0.45));

    % Baseline early to late
    n_baseline_licks_omitfirst = cell2mat(cellfun(@numel, baseline_licks_omitfirst,'UniformOutput',0)).*(1./0.45);
    lick_params.baseline_rate_mean_omitfirst_early = mean(n_baseline_licks_omitfirst(1:1000));
    lick_params.baseline_rate_mean_omitfirst_late = mean(n_baseline_licks_omitfirst(end-1000+1:end));
    lick_params.baseline_rate_shaping = (lick_params.baseline_rate_mean_omitfirst_early-lick_params.baseline_rate_mean_omitfirst_late)./(lick_params.baseline_rate_mean_omitfirst_early+lick_params.baseline_rate_mean_omitfirst_late);

end
