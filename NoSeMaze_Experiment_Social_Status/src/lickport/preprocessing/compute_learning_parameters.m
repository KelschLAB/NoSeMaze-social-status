function lick_params = compute_learning_parameters(d,cohort,plotPath)

% Loop over each animal's data
for an = 1:numel(d)

    % Pre-clearing
    clear lick_params_temp

    % Extract lickport data for this animal
    lickport_data = d(an).events;

    % Store the animal ID and cohort for later use
    lick_params_temp.ID = lickport_data(1).ID;
    lick_params_temp.cohort = cohort;

    % Display
    disp(['Processing animal nr. ' num2str(an)]);

    % Perform quality check: All trials must have the same ID
    assert(isscalar(unique({lickport_data.ID})), 'All trials must have the same ID.');

    %% Remove data after phase cutoff
    % lickport_data = remove_data_after_phase_cutoff(lickport_data, phase_cutoff);

    %% Compute baseline lick rate (Hz)
    lick_params_temp = compute_baseline_rate(lickport_data,lick_params_temp);

    %% Compute impulsivity (Hz)
    lick_params_temp = compute_impulsivity(lickport_data,lick_params_temp);    

    %% Compute reward prediction
    lick_params_temp = compute_reward_prediction(lickport_data,lick_params_temp,plotPath);

    %% CS+ and CS- lick modulation during the last 150 trials before reversal
    % Note that we only investigate the first second of the odor period
    % to assess modulation.
    lick_params_temp = compute_CS_lick_modulation(lickport_data, lick_params_temp);

    %% CS+ ramping
    % This function calculates the lick rate modulation (ramping) during the CS+ (conditioned stimulus)
    % period before the reversal, specifically comparing licks during the first and second segments of the odor presentation (= ramping).
    lick_params_temp = compute_CSplus_ramping(lickport_data, lick_params_temp);

    %% CS+ detection speed
    % This function calculates the speed at which the animal starts detecting or responding to the CS+ odor
    % based on the lick data recorded after odor onset. We're interested in the first time the licking rate exceeds
    % `baseline_rate_mean + 0.5 Hz` (the threshold). This is judged
    % within the 1 second of odor presentation.
    lick_params_temp = compute_CSplus_detection_speed(lickport_data, lick_params_temp);

    %% CS- detection speed
    % This function calculates the speed at which the animal starts detecting or responding to the CS+ odor
    % based on the lick data recorded after odor onset. We're interested in the first time the licking rate exceeds
    % `baseline_rate_mean - 0.5 Hz` (the threshold). This is judged
    % within the 1 second of odor presentation.
    lick_params_temp = compute_CSminus_detection_speed(lickport_data, lick_params_temp);

    %% CS+ valuation
    % This function calculates the peak to baseline within the time window 0.7-1.3 ms, in which the animal is evaluating the 
    % odor and potentially adapting its licks
    lick_params_temp = compute_CSplus_valuation(lickport_data, lick_params_temp);

    %% Giving Up at US or CS
    % This function checks for "giving-up" episodes in animal behavior
    % during different phases of a conditioning experiment, focusing on
    % the trials before and after a reversal. The function calculates
    % whether the animal's licking behavior significantly decreases during
    % the reversal, indicating a "giving-up" response.
    lick_params_temp = compute_giving_up_episodes(lickport_data, lick_params_temp);

    %% Switching Flexibility
    % Quantifies CS re-learning by calculating latency to switch (in number of trials)
    % This function computes the latency to switch behavior after a reversal phase,
    % by quantifying how quickly the animal adapts its response to the CS+ (rewarded)
    % and CS- (non-rewarded) stimuli. It uses the lick data before and after the reversal phase.
    % The thresholds are (for CS: 0.5-2.5s, US: 2.5-3.5s):
    % CS+ --> trials until reaching >70% of pre-reversal CS+ lick rate in
    % 6 out of 10 consecutive trials (moving average)
    % CS- --> trials until reducing to <50% of pre-reversal CS+ lick ratein
    % 6 out of 10 consecutive trials (moving average)
    % Also, switching latency shaping assesses the relation between switch latencies
    % in the first and the last reversal period.
    lick_params_temp = compute_switching_latency(lickport_data, lick_params_temp);

    %% Delay avoidance learner
    % This function calculates whether the animal has learned to avoid licking
    % during CS- trials (non-rewarded trials), which is interpreted as the
    % animal successfully learning the delay avoidance task.
    % It checks if the animal has made fewer than or equal to 1 lick, on average,
    % across all CS- trials during the CS period after odor delivery (0.5-2.5s).
    % Also, delay avooidance shaping assesses the relation between delay
    % avoidance in the first and the last reversal period.
    lick_params_temp = compute_delay_avoidance(lickport_data, lick_params_temp);

    %% Compute_pause_duration
    % This function identifies "giving-up" episodes where the animal shows a
    % drastic reduction in licking behavior during CS+ or US trials after a reversal.
    % The licking behavior is assessed by calculating lick counts before and after
    % a reversal phase for CS+ (rewarded) trials.
    lick_params_temp = compute_pause_duration(lickport_data, lick_params_temp);

    %% Correct hit and rejection rate in the 150 trials of CS+ and CS- before reversal and overall
    lick_params_temp = compute_hit_and_rejection_rates(lickport_data, lick_params_temp);

    %% Compute time to criterion of 80% is reached
    lick_params_temp = compute_time_to_criterion(lickport_data, lick_params_temp);

    %%
    % Find missing fields in lick_params
    if exist('lick_params', 'var')
        missing_fields = setdiff(fieldnames(lick_params), fieldnames(lick_params_temp));

        % Add missing fields to lick_params with default NaN values
        for i = 1:length(missing_fields)
            [lick_params_temp(:).(missing_fields{i})] = deal(NaN);  % Assign NaN to all existing elements
        end

        missing_fields = setdiff(fieldnames(lick_params_temp), fieldnames(lick_params));

        % Add missing fields to lick_params with default NaN values
        for i = 1:length(missing_fields)
            [lick_params(:).(missing_fields{i})] = deal(NaN);  % Assign NaN to all existing elements
        end
    end

    % Now assign lick_params_temp safely
    lick_params(an) = lick_params_temp;
end
end


