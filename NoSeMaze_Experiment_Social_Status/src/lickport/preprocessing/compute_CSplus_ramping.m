function lick_params = compute_CSplus_ramping(lickport_data, lick_params)
%% Compute CS+ Lick Rate Ramping
% This function calculates the ramping of the CS+ anticipation period.
% It compares the lick rates in two different time windows: the first
% (0.5s - 1.5s after odor onset) and the second (1.5s - 2.5s after odor onset).

% Extract the CS+ trials before reversal from the lick_params structure
cs_plus_trials_before_reversal = lick_params.cs_plus_trials_before_reversal;

% Initialize an empty array to store lick counts for each trial
cs_plus_lick_count = [];

% Define the time window for the odor presentation period (in seconds)
edges = [0.6,1.5,2.5];  % Licking counts from 0.6-1.5s and 1.5-2.5s

% Loop through each CS+ trial and calculate the lick counts
for tr = 1:numel(cs_plus_trials_before_reversal)
    % Histcounts counts the number of licks in the defined edges for each trial
    cs_plus_lick_count = cat(1, cs_plus_lick_count, ...
        histcounts(lickport_data(cs_plus_trials_before_reversal(tr)).licks_aft_od, edges));
end

% Compute the lick rate ramping by comparing the mean lick count in the second window to the first window
lick_params.cs_plus_ramping = mean(cs_plus_lick_count(:, 2)) / mean(cs_plus_lick_count(:, 1).*1/0.9);
end
