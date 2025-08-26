function [hierarchy_data] = extract_tube_competitions_from_LOG_clean(inputFile,double_detection_threshold)
%% extracts the tube test encounters from the RFID logs

% load input-file: LOG data containing detections of RFIDs
load(inputFile,'Data'); % loaded file 'Data' contains timestamps for every detection

% Only use the data from the detectors at the tube, not at the lickport (==0)
Data = Data(Data.Detector~=0,:);

ID = unique(Data.Animal);
ID(cellfun(@isempty, ID)) = [];
event_counter = 1;

%% remove double detections: on the same detector and less than double_detection_threshold (e.g, 0.5 seconds) from next detection
% this is a reasonably low threshold, maybe still some double detections
% in the Data but those can be dealt with more easily later now

Data.raw_index = [1:size(Data,1)]';
% add tube information
Data.tube = ones(size(Data,1),1);
Data.tube(Data.Detector==3 | Data.Detector==4) = 2;

Data = sortrows(Data,'Animal','ascend'); % sort Data by animals for double detection
delta_detection = diff(Data.Time_sec);
double_detections = [diff(Data.Detector)==0 & delta_detection<double_detection_threshold & strcmp(Data.Animal(1:end-1),Data.Animal(2:end)); false];
Data_filtered = Data(~double_detections,:);

% start_index = [true;~double_detections(1:end-1)];
start_index = find([true;~double_detections]);
Data_filtered.Detection_start = Data.Time_sec(start_index(1:size(Data_filtered,1)));

Data_filtered = sortrows(Data_filtered,'Time_sec','ascend'); % resort by time
Data = sortrows(Data,'Time_sec','ascend'); % resort by time


figure('visible','off');
tiledlayout('flow')
for animals = 1:length(ID)
    nexttile;
    data_this_animal = Data(contains(Data.Animal,ID{animals}),:);
    time_diff = diff(data_this_animal.Time_sec);
    histogram(time_diff(time_diff<3),'BinWidth',0.1,'Normalization','probability')
    line([.5 .5],[0 1],'LineStyle','--')
    xlim([0 3])
    title(ID{animals})
end
%% Loop over animal pairings (winners in rows, losers in columns)
% initialize output
match_matrix = zeros(length(ID),length(ID)); 
match_info = cell(length(ID),length(ID)); 
for winner = 1:length(ID)
    % find episodes of winner in the tube: assuming winner passes
    % through tube. If winner is in tube and defeats one or more opponents
    % but retracts at some point, difficult to capture. Think about making an exception case for
    % this.2
    winner_detections = Data_filtered(contains(Data_filtered.Animal,ID{winner}),:);
    
    % add original index information
    winner_detections.filtered_index = find(contains(Data_filtered{:,5},ID{winner}));

    % find winner passing through tube
    winner_in_tube = find(abs(diff(winner_detections.Detector))==1 & abs(diff(winner_detections.tube))==0); % JR: 29.04.2025
    
    % note: it is possible to look at tube preference (12 or 34
    % tube). Interesting?
    for events = 1:numel(winner_in_tube)

        winner_entry_time = winner_detections.Detection_start(winner_in_tube(events));
        winner_exit_time = winner_detections.Time_sec(winner_in_tube(events)+1);
        winner_time_spent_in = winner_exit_time-winner_entry_time; % in seconds
        winner_entry_Detector = winner_detections.Detector(winner_in_tube(events));
        winner_exit_Detector = winner_detections.Detector(winner_in_tube(events)+1);
        winner_tube = winner_detections.tube(winner_in_tube(events));
        multiple_animals_flag = false;
        
        % time threshold: 60 seconds arbitrary to exclude misdetections
        if winner_time_spent_in > 60; continue; end
        
        % loop over (putative) loser animals and check if animal was present
        % in the same tube while the (putative) winner was in
        for loser = 1:length(ID)
            if winner==loser; continue; end
            
            loser_detections = Data_filtered(contains(Data_filtered{:,5},ID{loser}),:);
            
            % critical step: find detections of loser while winner is in
            % tube AND at opposite entrance of winner AND loser also retracting
            % through that Detector again before the winner
            loser_exit_index = find(loser_detections.Time_sec>winner_entry_time ... 
                & loser_detections.Time_sec<winner_exit_time ...
                & winner_exit_Detector==loser_detections.Detector);
            
            if loser_exit_index % Note: JR (09.01.2025) changed it, as it might happen that the first detection (with loser_exit_index = 1) is already an event (e.g., Sarah's group 8, 2020-12-12)
                
                % exclude cases where winner is following loser
                if any(loser_exit_index > 1) % Note: JR (09.01.2025) changed it
                    if abs(loser_detections.Detector(loser_exit_index(end)-1)-loser_detections.Detector(loser_exit_index(end)))==1
                        continue;
                    end
                end
                
                % exclude event if time from start to end of defeat > 15
                % seconds
                if winner_time_spent_in>15
                    continue; 
                end
                
                % exclude event if more than two animals were involved
                if numel(unique(Data.Animal(find(Data.Time_sec > winner_entry_time & ...
                    Data.Time_sec < winner_exit_time & Data.tube==winner_tube))))>2
                    disp('multiple animals');
                    multiple_animals_flag = true;
%                     continue; 
                end
                
                % write 
                match_matrix(winner,loser) = match_matrix(winner,loser)+1;
                if isempty(match_info{winner,loser}); match_info{winner,loser}=struct; end
                dyad_counter =  match_matrix(winner,loser);
                match_info{winner,loser}(dyad_counter).winner_entry_time = winner_entry_time;
                match_info{winner,loser}(dyad_counter).winner_exit_time = winner_exit_time;
                match_info{winner,loser}(dyad_counter).winner_entry_Detector = winner_entry_Detector;
                match_info{winner,loser}(dyad_counter).winner_exit_Detector = winner_exit_Detector;
                match_info{winner,loser}(dyad_counter).winner_time_spent_in = winner_time_spent_in;
                match_info{winner,loser}(dyad_counter).winner_ID = ID{winner};
                match_info{winner,loser}(dyad_counter).multiple_animals = multiple_animals_flag;
                match_info{winner,loser}(dyad_counter).unique_event_ID = str2double([erase(Data.Day{1},{'-',' '}),num2str(event_counter)]);
                
                match_info{winner,loser}(dyad_counter).loser_exit_time = loser_detections.Time_sec(loser_exit_index);
                match_info{winner,loser}(dyad_counter).loser_exit_Detector = loser_detections.Detector(loser_exit_index);
                match_info{winner,loser}(dyad_counter).loser_time_detected = loser_detections.Time_sec(loser_exit_index)-loser_detections.Detection_start(loser_exit_index);
                disp(['W: ' ID{winner} '; L: ' ID{loser}])
                if length(loser_detections.Detection_start(loser_exit_index))>1
                    match_info{winner,loser}(dyad_counter).loser_time_spent_in = loser_detections.Time_sec(loser_exit_index(end)) - loser_detections.Detection_start(loser_exit_index(1));
                % elseif loser_exit_index>1 && loser_detections.Detection_start(loser_exit_index-1) > (loser_detections.Time_sec(loser_exit_index)-10)
                elseif loser_exit_index>1 && loser_detections.Detection_start(loser_exit_index-1) > (loser_detections.Time_sec(loser_exit_index)-10) && loser_detections.Detector(loser_exit_index)~=loser_detections.Detector(loser_exit_index-1)
                    warning('detectors DIFFERENT');
                    match_info{winner,loser}(dyad_counter).loser_time_spent_in = loser_detections.Time_sec(loser_exit_index) - loser_detections.Detection_start(loser_exit_index);
                elseif loser_exit_index>1 && loser_detections.Detection_start(loser_exit_index-1) > (loser_detections.Time_sec(loser_exit_index)-10) && loser_detections.Detector(loser_exit_index)==loser_detections.Detector(loser_exit_index-1) % REPLACED: 19.05.2025, JR
                    match_info{winner,loser}(dyad_counter).loser_time_spent_in = loser_detections.Time_sec(loser_exit_index) - loser_detections.Detection_start(loser_exit_index-1);
                else
                    match_info{winner,loser}(dyad_counter).loser_time_spent_in = loser_detections.Time_sec(loser_exit_index) - loser_detections.Detection_start(loser_exit_index);
                end
                match_info{winner,loser}(dyad_counter).loser_ID = ID{loser};

                % add filtered data snippet
                match_info{winner,loser}(dyad_counter).data_filtered = Data_filtered(max([1,winner_detections.filtered_index(winner_in_tube(events))-20]):...
                    min([size(Data_filtered,1),winner_detections.filtered_index(winner_in_tube(events))+20]),:);
                
                % add raw data snippet
                match_info{winner,loser}(dyad_counter).data = Data(max([1,winner_detections.raw_index(winner_in_tube(events))-20]):...
                    min([size(Data_filtered,1),winner_detections.raw_index(winner_in_tube(events))+20]),:);
                
                event_counter = event_counter+1;
                
            end
        end    
    end
end


%% parse output
hierarchy_data.ID = ID;
hierarchy_data.match_matrix = match_matrix;
hierarchy_data.match_info = match_info;
hierarchy_data.Data = Data;
hierarchy_data.Data_filtered = Data_filtered;
hierarchy_data.input = inputFile;
end

