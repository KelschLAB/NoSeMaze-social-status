function [d] = clean_phases_and_exclude_animals(d,cohort,exclusion_table)


% Exclude animals
for an=1:length(d)
    if any(strcmp(exclusion_table.cohort, cohort) & strcmp(exclusion_table.animalID, d(an).events(1).ID))
        exclusion_vector(an)=true;
    else
        exclusion_vector(an)=false;
    end
end
% Exclusion
d=d(~exclusion_vector);

% loop over animals
for an=1:length(d)

    % Read phases for current animal
    current_phases = [d(an).events.phase];

    % Use tabulate to count the occurrences
    current_table = tabulate(current_phases);

    % Display information
    disp([cohort ', ' d(an).events(1).ID ': ' num2str(sum(current_table(:,2)<50)) ' phases with less then 50 trials' ])

    % Define phases to delete
    phases_to_delete = current_table(current_table(:,2)<50,1);
    % Delete phases
    d(an).events = d(an).events(~ismember([d(an).events.phase],phases_to_delete));
    % Display
    disp([num2str(length(phases_to_delete)) ' phases deleted.']);

    % Rename phases
    phases_to_keep = current_table(current_table(:,2)>=50,1);
    renaming_phases_table = table(phases_to_keep,[1:length(phases_to_keep)]','VariableNames',{'OldPhaseName','NewPhaseName'});
    for phase_id = 1:height(renaming_phases_table)
        % Find the indices of the events where phase matches OldPhaseName
        indices_to_update = [d(an).events.phase] == renaming_phases_table.OldPhaseName(phase_id);
        % Update the phase for those specific events
        [d(an).events(indices_to_update).phase] = deal(renaming_phases_table.NewPhaseName(phase_id));
    end
end