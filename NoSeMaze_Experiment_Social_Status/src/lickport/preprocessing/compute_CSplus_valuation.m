function lick_params = compute_CSplus_valuation(lickport_data, lick_params)
    %% CS+ valuation
    % This function calculates the speed at which the animal starts detecting or responding to the CS+ odor 
    % based on the lick data recorded after odor onset.

    % Extract the CS+ trials before reversal from the lick_params structure
    cs_plus_trials_before_reversal = lick_params.cs_plus_trials_before_reversal;

    % Initialize an empty array to store lick counts for each trial
    cs_plus_lick_count = [];

    % Define the time window of interest (from 0.7s to 2.5s after odor onset)
    % This time window will be used to bin the lick timestamps for analysis.
    % Start at 700 ms is used to a certain time lag after odor onset that
    % is expected to be needed for "valuation"
    edges = 0.7:0.05:1.3;  % Creates edges: [0.7, 0.75, 0.8, ..., 2.5]

    % Loop over each CS+ trial before reversal
    for tr = 1:numel(cs_plus_trials_before_reversal)
        % For each trial, calculate the high-temporal resolution PSTH (peristimulus time histogram)
        % This counts the number of licks that occur within the specified time bins defined by `edges`.
        cs_plus_lick_count = cat(1, cs_plus_lick_count, ...
            histcounts(lickport_data(cs_plus_trials_before_reversal(tr)).licks_aft_od, edges));
    end

    % Compute the average lick rate (PSTH) over the trials
    % We multiply the result by 20 to adjust the temporal resolution.
    % The reason for omitting the first bin is to exclude the initial lick burst typically associated with odor delivery.
    psth = mean(cs_plus_lick_count) * 20;  % This is the average lick rate over the time window.

    % CS+ detection speed: 
    % This calculates the time (in seconds) at which the licking rate exceeds the baseline by a certain threshold.
    % We're interested in the first time the licking rate exceeds `baseline_rate_mean + 0.5 Hz` (the threshold).
    % The `find` function locates the first bin that exceeds this threshold.
    [max_val,max_ind] = max(psth);
    lick_params.cs_plus_valuation_peak = max_val;
    lick_params.cs_plus_valuation_peak_to_base = max_val./lick_params.baseline_rate_CSplus_mean_omitfirst;
    lick_params.cs_plus_valuation_peak_minus_base = max_val-lick_params.baseline_rate_CSplus_mean_omitfirst;
    lick_params.cs_plus_valuation_time_to_peak = 0.2 + 0.05 * max_ind;

    % If the detection speed is not found (i.e., no threshold exceeded), set it to NaN
    if isempty(lick_params.cs_plus_detection_speed)
        lick_params.cs_plus_valuation_peak = NaN;
        lick_params.cs_plus_valuation_peak_to_base = NaN;
        lick_params.cs_plus_valuation_time_to_peak = NaN;
    end
end
