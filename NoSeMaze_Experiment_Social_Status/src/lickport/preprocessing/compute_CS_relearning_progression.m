function lick_params = compute_CS_relearning_progression(lick_params,phase_cutoff)
    %% Compute CS+ and CS- Relearning Progression
    % This function calculates the relearning progression for both CS+ and CS- trials
    % by computing the ratio of the switch latencies during the second reversal (phase 2)
    % to those from the last reversal (determined by phase_cutoff).
    %
    % The relearning progression is calculated as:
    %    (switch latency at reversal 2) / (switch latency at last reversal)
    %
    try
        % Calculate relearning progression
        lick_params.cs_plus_relearning_progression_at_cs = lick_params.cs_plus_switch_latency_at_cs_rev2/lick_params.(['cs_plus_switch_latency_at_cs_rev',num2str(phase_cutoff-1)]);
        lick_params.cs_plus_relearning_progression_at_us = lick_params.cs_plus_switch_latency_at_us_rev2/lick_params.(['cs_plus_switch_latency_at_us_rev',num2str(phase_cutoff-1)]);
        lick_params.cs_minus_relearning_progression_at_cs = lick_params.cs_minus_switch_latency_at_cs_rev2/lick_params.(['cs_minus_switch_latency_at_cs_rev',num2str(phase_cutoff-1)]);
        lick_params.cs_minus_relearning_progression_at_us = lick_params.cs_minus_switch_latency_at_us_rev2/lick_params.(['cs_minus_switch_latency_at_us_rev',num2str(phase_cutoff-1)]);
    catch
        % If any required data is missing or calculation fails, set to NaN
        lick_params.cs_plus_relearning_progression_at_cs = NaN;
        lick_params.cs_plus_relearning_progression_at_us = NaN;
        lick_params.cs_minus_relearning_progression_at_cs = NaN;
        lick_params.cs_minus_relearning_progression_at_us = NaN;
    end
end
