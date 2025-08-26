function lick_params = compute_impulsivity(lickport_data,lick_params)
    % This function computes the lick rate during the first 100 ms
    % (at trial, odor, reward)

    % Case 1: Average lick-rate (Hz) during baseline from 0.0 to 0.1s after trial-start
    lick_peak_trial = cellfun(@(x) x(x<0.1),{lickport_data(:).licks_bef_od},'UniformOutput',0);
    lick_params.impulsivity_trialStart = mean(cell2mat(cellfun(@numel, lick_peak_trial,'UniformOutput',0)).*(1./0.1));
    lick_params.impulsivity_trialStart_norm = (lick_params.impulsivity_trialStart-lick_params.baseline_rate_mean_omitfirst)./(lick_params.impulsivity_trialStart+lick_params.baseline_rate_mean_omitfirst);

    % Case 2: Average lick-rate (Hz) during odor from 0.5 to 0.6s after trial-start
    lick_peak_odor = cellfun(@(x) x(x>0.5 & x<0.6),{lickport_data(:).licks_aft_od},'UniformOutput',0);
    lick_params.impulsivity_odorStart = mean(cell2mat(cellfun(@numel, lick_peak_odor,'UniformOutput',0)).*(1./0.1));
    lick_params.impulsivity_odorStart_norm = (lick_params.impulsivity_odorStart-lick_params.baseline_rate_mean_omitfirst)./(lick_params.impulsivity_odorStart+lick_params.baseline_rate_mean_omitfirst);

    % Case 3: Average lick-rate (Hz) during reward from 2.5 to 2.6 s after trial-start
    lick_peak_reward = cellfun(@(x) x(x>2.5 & x<2.6),{lickport_data(:).licks_aft_od},'UniformOutput',0);
    lick_params.impulsivity_rewardStart = mean(cell2mat(cellfun(@numel, lick_peak_reward,'UniformOutput',0)).*(1./0.1));
    lick_params.impulsivity_rewardStart_norm = (lick_params.impulsivity_rewardStart-lick_params.baseline_rate_mean_omitfirst)./(lick_params.impulsivity_rewardStart+lick_params.baseline_rate_mean_omitfirst);
end
