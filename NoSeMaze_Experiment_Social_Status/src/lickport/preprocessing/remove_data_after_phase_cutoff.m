function lickport_data = remove_data_after_phase_cutoff(lickport_data, phase_cutoff)
    %% Remove data after phase cutoff
    % This function removes trials where the phase exceeds the given phase_cutoff.
    % The cutoff is used to filter out data after a certain phase of the experiment.
    
    % Remove trials where the phase exceeds phase_cutoff
    lickport_data([lickport_data.phase]' > phase_cutoff) = [];
end
