function lick_params = compute_CSplus_detection_speed(lickport_data, lick_params)
    %% CS+ detection/valuation speed
    % This function calculates the speed at which the animal starts detecting or responding to the CS+ odor 
    % based on the lick data recorded after odor onset.

    % Extract the CS+ trials before reversal from the lick_params structure
    cs_plus_trials_before_reversal = lick_params.cs_plus_trials_before_reversal;

    % Initialize an empty array to store lick counts for each trial
    cs_plus_lick_count = [];

    % Define the time window of interest (from 0.7s to 2.5s after odor onset)
    % This time window will be used to bin the lick timestamps for analysis.
    edges = 0.5:0.05:2.5;  % Creates edges: [0.7, 0.75, 0.8, ..., 2.5]

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
    % The first 100 ms are not considered. 
    lick_params.cs_plus_detection_speed = 0.1 + 0.05 * find(psth(3:end) > 0.5 + lick_params.baseline_rate_CSplus_mean_omitfirst, 1, 'first');

    % If the detection speed is not found (i.e., no threshold exceeded), set it to NaN
    if isempty(lick_params.cs_plus_detection_speed)
        lick_params.cs_plus_detection_speed = NaN;
    end
end
