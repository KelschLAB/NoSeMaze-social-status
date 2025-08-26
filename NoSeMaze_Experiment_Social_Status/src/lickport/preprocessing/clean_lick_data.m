function d = clean_lick_data(d)
%% Data Cleaning: minimum lick-timestamp difference = 50ms
% This function checks for consecutive licks within 50ms and removes them.
% It ensures that only valid licks (50ms apart) are kept.

% Loop over each animal's data
for an = 1:numel(d)
    % Extract lickport data for this animal
    lickport_data = d(an).events;

    % Iterate over each trial to clean the lick data
    for tr = 1:numel(lickport_data)
        % Concatenate before and after odor licks into one array
        all_trial_licks = cat(2, lickport_data(tr).licks_bef_od, lickport_data(tr).licks_aft_od);

        % Skip the trial if there are no licks recorded
        if isempty(all_trial_licks)
            disp(['Trial ' num2str(tr) ' has no lick data. Skipping.']);
            continue;  % Skip this trial and move to the next one
        end

        % Sort the licks in ascending order
        all_trial_licks = sort(all_trial_licks);

        % Initialize valid licks with the first lick in the trial
        valid_licks = all_trial_licks(1);

        % Loop over the rest of the licks to remove any that are less than 50ms apart
        for lc = 2:numel(all_trial_licks)
            if all_trial_licks(lc) - valid_licks(end) >= 0.05  % If the difference is >= 50ms
                valid_licks = [valid_licks, all_trial_licks(lc)];
            end
        end

        if length(all_trial_licks)~=length(valid_licks)
            disp([num2str(length(all_trial_licks)-length(valid_licks)) ' lick(s) deleted']);
        end

        % Split licks back into before and after odor
        lickport_data(tr).licks_bef_od = valid_licks(valid_licks < 0.5);
        lickport_data(tr).licks_aft_od = valid_licks(valid_licks >= 0.5);
    end

    d(an).events = lickport_data;
end
end
