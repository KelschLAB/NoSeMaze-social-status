function full_hierarchy = combine_days_for_hierarchy(inputPath,interimOutputPath,cohort,start_date,end_date)
% Load exclude, invert, and include events
exclude_file_path = fullfile(inputPath,'exclude_events.mat');
if isfile(exclude_file_path); load(exclude_file_path,'exclude'); else; exclude = 0; end
if size(exclude,1)<=size(exclude,2); exclude=exclude'; end % assert orientation of arrray
invert_file_path = fullfile(inputPath,'invert_events.mat');
if isfile(invert_file_path); load(invert_file_path,'invert'); else; invert = 0; end
if size(invert,1)<=size(invert,2); invert=invert'; end % assert orientation of arrray
include_file_path = fullfile(inputPath,'include_events.mat');
if isfile(include_file_path); load(include_file_path,'include'); else; include = 0; end
if size(include,1)<=size(include,2); include=include'; end % assert orientation of arrray

% Get all expected days from the cohort's start and end dates
all_days = datestr((datetime(start_date):datetime(end_date))', 'yyyy-mm-dd');

% Find existing Data.mat files
Data_list = dir(fullfile(inputPath,'*',filesep,'Data.mat'));
existing_days = arrayfun(@(x) x.folder(end-9:end), Data_list, 'UniformOutput', false); % Extract dates from filenames
missing_days = setdiff(all_days, existing_days);

% Initialize full_hierarchy
full_hierarchy = struct();

% first day with existing hierarchy
reference_firstDataDay = find(strcmp(all_days,existing_days(1)));

% Loop over all expected days (existing and missing)
for day_num = 1:size(all_days,1)
    cur_day = all_days(day_num, :);
    
    % Check if current day has a Data.mat file
    data_file_idx = find(strcmp(existing_days, cur_day), 1);
    if isempty(data_file_idx)        
        % Create a placeholder hierarchy_data for missing days
        hierarchy_data.ID = {};
        hierarchy_data.match_matrix = [];
        hierarchy_data.match_info = {};
        hierarchy_data.Data = [];
        hierarchy_data.Data_filtered = [];
        hierarchy_data.input = [];
    elseif ~isempty(data_file_idx)   
        % Load hierarchy_data for the current day
        load(fullfile(Data_list(data_file_idx).folder, Data_list(data_file_idx).name), 'hierarchy_data');
        
        %%% to Do: account for decrease in animal number or swap of
        %%% RFID identitiy!
        cur_IDs = hierarchy_data.ID;
        if day_num>reference_firstDataDay

            % assume no new animals/RFIDs
            if any(contains(cur_IDs,full_hierarchy(reference_firstDataDay).ID)==0) % for some reason misdetection of IDs that do not belong to cohort -> exclude
                delete_animal_index = find(contains(cur_IDs,full_hierarchy(reference_firstDataDay).ID)==0);
                hierarchy_data.ID(delete_animal_index) = [];
                hierarchy_data.match_info(:,delete_animal_index) = [];
                hierarchy_data.match_info(delete_animal_index,:) = [];
                hierarchy_data.match_matrix(:,delete_animal_index) = [];
                hierarchy_data.match_matrix(delete_animal_index,:) = [];
                cur_IDs = hierarchy_data.ID;
            end
            
            assert(~any(contains(cur_IDs,full_hierarchy(reference_firstDataDay).ID)==0))
            
            
            % if animal was taken out, resort hierarchy_data so that the
            % matrices match the full_hierarchy structure
            if numel(cur_IDs)<numel(full_hierarchy(reference_firstDataDay).ID)
                % parse original data
                full_hierarchy(day_num).ID_misaligned = hierarchy_data.ID;
                full_hierarchy(day_num).match_matrix_misaligned = hierarchy_data.match_matrix;
                full_hierarchy(day_num).match_info_misaligned = hierarchy_data.match_info;
                
                % rebuild match_matrix and match_info in correct order
                match_matrix_realigned = zeros(numel(full_hierarchy(reference_firstDataDay).ID),numel(full_hierarchy(reference_firstDataDay).ID));
                match_info_realigned = cell(numel(full_hierarchy(reference_firstDataDay).ID),numel(full_hierarchy(reference_firstDataDay).ID));
                IDs_realigned = full_hierarchy(reference_firstDataDay).ID;
                IDs_realigned{~contains(full_hierarchy(reference_firstDataDay).ID, hierarchy_data.ID)} = [];
                
                % loop over misaligned data, find non-zero entries and parse to
                % correct position in matrix
                for i=1:size(hierarchy_data.match_matrix,1)
                    for j=1:size(hierarchy_data.match_matrix,2)
                        if hierarchy_data.match_matrix(i,j)>0
                            for n=1:size(hierarchy_data.match_info{i,j},2)
                                correct_winner_position = find(contains(full_hierarchy(reference_firstDataDay).ID,hierarchy_data.match_info{i,j}(n).winner_ID));
                                correct_loser_position = find(contains(full_hierarchy(reference_firstDataDay).ID,hierarchy_data.match_info{i,j}(n).loser_ID));
                                
                                % write to match_matrix
                                match_matrix_realigned(correct_winner_position,correct_loser_position)=match_matrix_realigned(correct_winner_position,correct_loser_position)+1;
                                
                                % write to match_info
                                match_info_realigned{correct_winner_position,correct_loser_position}(n)=hierarchy_data.match_info{i,j}(n);
                                
                            end
                        end
                    end
                end
                
                % parse realigned data to hierarchy_data
                hierarchy_data.match_matrix = match_matrix_realigned;
                hierarchy_data.match_info = match_info_realigned;
                hierarchy_data.ID = full_hierarchy(reference_firstDataDay).ID;
                hierarchy_data.ID = IDs_realigned;
            else % else assert same sorting
                assert(strcmp(cat(1,cur_IDs{:}),cat(1,full_hierarchy(reference_firstDataDay).ID{:})))
            end
        end
        
        
        % exclude multi-animal events (not if manually re-included or inverted)
        for i=1:size(hierarchy_data.match_matrix,1)
            for j=1:size(hierarchy_data.match_matrix,2)
                if ~isempty(hierarchy_data.match_info{i,j})
                    to_delete = [];
                    for n=1:size(hierarchy_data.match_info{i,j},2)
                        if hierarchy_data.match_info{i,j}(n).multiple_animals && ~ismember(hierarchy_data.match_info{i,j}(n).unique_event_ID, include) && ~ismember(hierarchy_data.match_info{i,j}(n).unique_event_ID, invert)
                            hierarchy_data.match_matrix(i,j)=hierarchy_data.match_matrix(i,j)-1;
                            to_delete = cat(1,to_delete,n);
                        end
                    end
                    hierarchy_data.match_info{i,j}(to_delete)=[];
                end
            end
        end
        
        % exclude manually-excluded events
        for i=1:size(hierarchy_data.match_matrix,1)
            for j=1:size(hierarchy_data.match_matrix,2)
                if ~isempty(hierarchy_data.match_info{i,j})
                    to_delete = [];
                    for n=1:size(hierarchy_data.match_info{i,j},2)
                        if ismember(hierarchy_data.match_info{i,j}(n).unique_event_ID, exclude)
                            hierarchy_data.match_matrix(i,j)=hierarchy_data.match_matrix(i,j)-1;
                            to_delete = cat(1,to_delete,n);
                        end
                    end
                    hierarchy_data.match_info{i,j}(to_delete)=[];
                end
            end
        end
        
        %%% From here changed by Jonathan Reinwald, 12.12.2024
        % %         % invert winner-loser relationship in manually-curated events (e.g.,
        % %         % because of back-and-forth pushing the automatic extraction
        % %         % incorrectly captures the events)
        % %         for i=1:size(hierarchy_data.match_matrix,1)
        % %             for j=1:size(hierarchy_data.match_matrix,2)
        % %                 if ~isempty(hierarchy_data.match_info{i,j})
        % %                     to_delete = [];
        % %                     for n=1:size(hierarchy_data.match_info{i,j},2)
        % %                         if ismember(hierarchy_data.match_info{i,j}(n).unique_event_ID, invert)
        % %                             hierarchy_data.match_matrix(i,j)=hierarchy_data.match_matrix(i,j)-1;
        % %                             hierarchy_data.match_matrix(j,i)=hierarchy_data.match_matrix(j,i)+1;
        % %                             hierarchy_data.match_info{j,i}(numel(hierarchy_data.match_info{j,i})+1)=hierarchy_data.match_info{i,j}(n);
        % %                             to_delete = cat(1,to_delete,n);
        % %                         end
        % %                     end
        % %                     hierarchy_data.match_info{i,j}(to_delete)=[];
        % %                 end
        % %             end
        % %         end
        
        % Ensure all match_info entries have a field to track inversion
        for i = 1:size(hierarchy_data.match_info, 1)
            for j = 1:size(hierarchy_data.match_info, 2)
                if ~isempty(hierarchy_data.match_info{i, j})
                    for n = 1:size(hierarchy_data.match_info{i, j}, 2)
                        hierarchy_data.match_info{i, j}(n).already_inverted = false; % Initialize as not inverted
                    end
                end
            end
        end
        
        
        
        % Loop through match_matrix to invert winner-loser relationships
        for i = 1:size(hierarchy_data.match_matrix, 1)
            for j = 1:size(hierarchy_data.match_matrix, 2)
                if ~isempty(hierarchy_data.match_info{i, j})
                    to_delete = [];
                    for n = 1:size(hierarchy_data.match_info{i, j}, 2)
                        % Check if the event matches the invert list and hasn't been inverted yet
                        if ismember(hierarchy_data.match_info{i, j}(n).unique_event_ID, invert) && ...
                                ~hierarchy_data.match_info{i, j}(n).already_inverted
                            
                            % Update match_matrix to reflect the inversion of the match result
                            hierarchy_data.match_matrix(i, j) = hierarchy_data.match_matrix(i, j) - 1;
                            hierarchy_data.match_matrix(j, i) = hierarchy_data.match_matrix(j, i) + 1;
                            
                            % Flip winner and loser IDs
                            new_loser = hierarchy_data.match_info{i, j}(n).winner_ID;
                            new_winner = hierarchy_data.match_info{i, j}(n).loser_ID;
                            hierarchy_data.match_info{i, j}(n).winner_ID = new_winner;
                            hierarchy_data.match_info{i, j}(n).loser_ID  = new_loser;
                            
                            % -------------------------------------------------
                            % ADJUST TIME-RELATED DETAILS FOR THE INVERTED MATCH
                            % -------------------------------------------------
                            
                            % Locate the row in data_filtered where the old loser's exit time is recorded
                            index_old_loser_exit = find(abs( ...
                                hierarchy_data.match_info{i,j}(n).data_filtered.Time_sec - ...
                                hierarchy_data.match_info{i,j}(n).loser_exit_time(end)) < 1e-6);
                            
                            % Save details of the old loser for updating new winner's entry details
                            old_loser_exit_time      = hierarchy_data.match_info{i,j}(n).loser_exit_time(end);
                            old_loser_exit_Detector  = hierarchy_data.match_info{i,j}(n).loser_exit_Detector;
                            
                            % Save details of the old winner for updating new loser's exit details
                            old_winner_exit_Detector  = hierarchy_data.match_info{i,j}(n).winner_exit_Detector;
                            old_winner_entry_Detector = hierarchy_data.match_info{i,j}(n).winner_entry_Detector;
                            
                            % Update the new winner's entry time to the old loser's exit time
                            hierarchy_data.match_info{i,j}(n).winner_entry_time = old_loser_exit_time;
                            
                            % Find indices in data_filtered where the new winner is detected
                            % at the old winner's entry detector
                            indices_new_winner_at_new_exit = find( ...
                                strcmp(hierarchy_data.match_info{i,j}(n).data_filtered.Animal, new_winner) & ...
                                (hierarchy_data.match_info{i,j}(n).data_filtered.Detector == old_winner_entry_Detector));
                            
                            % Filter for rows where the detection occurs after the old loser's exit time
                            indices_new_winner_at_new_exit = ...
                                indices_new_winner_at_new_exit(indices_new_winner_at_new_exit > index_old_loser_exit);
                            
                            % Assign the new winner's entry detector to the old winner's exit detector
                            hierarchy_data.match_info{i,j}(n).winner_entry_Detector = old_winner_exit_Detector;
                            
                            % Assign the new winner's exit time using the first valid detection found
                            if ~isempty(indices_new_winner_at_new_exit)
                                hierarchy_data.match_info{i,j}(n).winner_exit_time = ...
                                    hierarchy_data.match_info{i,j}(n).data_filtered.Time_sec(indices_new_winner_at_new_exit(1));
                            else
                                % If no new winner exit is found, fallback logic or skip
                                continue;
                            end
                            
                            % Assign the new winner's exit detector to the old winner's entry detector
                            hierarchy_data.match_info{i,j}(n).winner_exit_Detector = old_winner_entry_Detector;
                            
                            % Calculate the time the new winner spent in the detector
                            hierarchy_data.match_info{i,j}(n).winner_time_spent_in = ...
                                hierarchy_data.match_info{i,j}(n).winner_exit_time - ...
                                hierarchy_data.match_info{i,j}(n).winner_entry_time;
                            
                            % Find indices in data_filtered where the new loser is detected at the old winner's entry detector
                            indices_new_loser_at_new_exit = find( ...
                                strcmp(hierarchy_data.match_info{i,j}(n).data_filtered.Animal, new_loser) & ...
                                (hierarchy_data.match_info{i,j}(n).data_filtered.Detector == old_winner_entry_Detector));
                            
                            % Filter for rows where the detection occurs before the new winner's exit time
                            indices_new_loser_at_new_exit = ...
                                indices_new_loser_at_new_exit(indices_new_loser_at_new_exit < indices_new_winner_at_new_exit(1));
                            
                            % Assign the new loser's exit time using the last valid detection found
                            if ~isempty(indices_new_loser_at_new_exit)
                                hierarchy_data.match_info{i,j}(n).loser_exit_time = ...
                                    hierarchy_data.match_info{i,j}(n).data_filtered.Time_sec(indices_new_loser_at_new_exit(end));
                            else
                                % If no new loser exit is found, fallback logic or skip
                                continue;
                            end
                            
                            % Assign the new loser's exit detector to the old winner's entry detector
                            hierarchy_data.match_info{i,j}(n).loser_exit_Detector = old_winner_entry_Detector;
                            
                            % Assign NaN to the loser's time detected (if not calculated elsewhere)
                            hierarchy_data.match_info{i,j}(n).loser_time_detected = NaN;
                            
                            % -------------------------------------------------
                            % MOVE THE INVERTED EVENT TO THE [j, i] CELL
                            % -------------------------------------------------
                            
                            if isempty(hierarchy_data.match_info{j, i})
                                hierarchy_data.match_info{j, i} = hierarchy_data.match_info{i, j}(n);
                            else
                                hierarchy_data.match_info{j, i}(end + 1) = hierarchy_data.match_info{i, j}(n);
                            end
                            
                            % Mark the event as inverted in both places
                            hierarchy_data.match_info{j, i}(end).already_inverted = true;
                            hierarchy_data.match_info{i, j}(n).already_inverted   = true;
                            
                            % Mark for deletion from the original position
                            to_delete = cat(1, to_delete, n);
                        end
                    end
                    % Delete inverted events from the original position
                    hierarchy_data.match_info{i, j}(to_delete) = [];
                end
            end
        end
    end
        
    % parse to combined data structure
    full_hierarchy(day_num).ID = hierarchy_data.ID;
    full_hierarchy(day_num).match_matrix = hierarchy_data.match_matrix;
    full_hierarchy(day_num).match_info = hierarchy_data.match_info;
    full_hierarchy(day_num).Data = hierarchy_data.Data;
    full_hierarchy(day_num).Data_filtered = hierarchy_data.Data_filtered;
    full_hierarchy(day_num).input = hierarchy_data.input;
    
    % assert all IDs are matching
    %     assert(mean(cell2mat(full_hierarchy(reference_firstDataDay).ID)==cell2mat(hierarchy_data.ID),'all')==1);
    
end

%% save output
save(fullfile(interimOutputPath,['full_hierarchy_' cohort '.mat']),'full_hierarchy');


