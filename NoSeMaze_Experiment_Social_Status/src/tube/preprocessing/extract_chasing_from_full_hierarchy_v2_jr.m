function full_hierarchy = extract_chasing_from_full_hierarchy_v2_jr(full_hierarchy,ops)
%% Extract chasing/close-following events from the tube detections
% David Wolf, 01/2023
%
%   Input: 
%       full_hierarchy
%       ops: specify options like thresholds, plotting functions etc

%% set defaults if not specified
if ~isfield(ops,'plotter')
    ops.plotter = 0;
end
if ~isfield(ops,'save_dir')
    ops.save_dir = cd;
end
if ~isfield(ops,'threshold_lag_at_detector')
    ops.threshold_lag_at_detector = 1.5; % in seconds
end
if ~isfield(ops,'threshold_lag_through_tube')
    ops.threshold_lag_through_tube = 3; % in seconds
end

%% refilter raw data to fit the data from 2023 better

for day = 1:size(full_hierarchy,2)
   
    % raw data from this day
    Data = full_hierarchy(day).Data;
    
    Data = sortrows(Data,'Animal','ascend'); % sort Data by animals for double detection
    delta_detection = diff(Data.Time_sec);
    double_detections = [diff(Data.Detector)==0 & delta_detection<1 & strcmp(Data.Animal(1:end-1),Data.Animal(2:end)); false];
    Data_filtered = Data(~double_detections,:);

    % start_index = [true;~double_detections(1:end-1)];
    start_index = find([true;~double_detections]);
    Data_filtered.Detection_start = Data.Time_sec(start_index(1:size(Data_filtered,1)));

    full_hierarchy(day).Data_filtered = sortrows(Data_filtered,'Time_sec','ascend'); % resort by time   
    
end


%%
% using the Data_filtered find events where one mouse is chasing the
% other (defined by the thresholds above)
filtered_detections = cat(1,full_hierarchy.Data_filtered);
counter = 1;
for day = 1:size(full_hierarchy,2)
    events_this_day = size(full_hierarchy(day).Data_filtered,1);
    filtered_detections.Detection_start(counter:counter+events_this_day-1) = filtered_detections.Detection_start(counter:counter+events_this_day-1)+(day-1)*24*60*60; % add offset for concating days
    counter = counter+events_this_day;
end
filtered_detections = sortrows(filtered_detections,'Detection_start');
filtered_detections(~contains(filtered_detections.Animal,full_hierarchy(1).ID),:) = [];

% find tube passings where animals are closely following each other:
% same detector % time_diff under threshold_dist 
% separate tubes from each other so the time diff doesnt interfere with
% each other
chasing_index = [];
for tt = 1:2
    tube_filtered = filtered_detections(filtered_detections.tube==tt,:);
    tube_filtered(tube_filtered.Detector==0,:)=[];


    % find putative chasing events: the sequence has to be 1-2-1-2
    % (animals) with certain lag criteria
    tmp1 = strfind(abs(diff(tube_filtered.Detector))',[0 1 0]);% abs is the difference to the former script; by JR 14.12.2023
    tmp2 = find(~strcmp(tube_filtered.Animal(1:end-3),tube_filtered.Animal(2:end-2)) & strcmp(tube_filtered.Animal(1:end-3),tube_filtered.Animal(3:end-1)) & ...
        ~strcmp(tube_filtered.Animal(3:end-1),tube_filtered.Animal(4:end)) & strcmp(tube_filtered.Animal(2:end-2),tube_filtered.Animal(4:end)))';
    close_following_index = tmp1(ismember(tmp1,tmp2));


    % time lags between mice at detectors
    delta_t_entry = tube_filtered.Detection_start(close_following_index+1)-tube_filtered.Detection_start(close_following_index);
    delta_t_exit = tube_filtered.Detection_start(close_following_index+3)-tube_filtered.Detection_start(close_following_index+2);

    % visualize temporal distribution of potential close following events
     if ops.plotter
        f=figure;
        s(1)=subplot(1,2,1);
        histogram(delta_t_entry(delta_t_entry<8),'BinWidth',.25)
        title('\Delta at entry')
        xlim([0 5]);
        ylabel('count'); xlabel('\Delta time chase to victim');
        set_fonts()
        
        s(2)=subplot(1,2,2);
        histogram(delta_t_exit(delta_t_exit<8),'BinWidth',.25)
        title('\Delta at exit')
        xlim([0 5]);
        ylabel('count'); xlabel('\Delta time chase to victim');
        linkaxes(s,'xy');
        set_fonts()
        
        f.Units = 'centimeters';
        f.Position = [3 3 8 4];
        exportgraphics(f,fullfile(ops.save_dir,['time_diff_at_detector',num2str(tt),'.png']));
    end
    
    % filter out following events above threshold for passing through
    % tube
    delta_through_victim = tube_filtered.Detection_start(close_following_index+2)-tube_filtered.Detection_start(close_following_index);
    delta_through_bully = tube_filtered.Detection_start(close_following_index+3)-tube_filtered.Detection_start(close_following_index+1);

    if ops.plotter
        f=figure;
        s(1)=subplot(1,2,1);
        histogram(delta_through_victim,'BinWidth',.1)
        title('\Delta through tube victim')
        xlabel('\Delta time through tube');
        ylabel('count');
        xlim([0 3]);
        set_fonts()
        
        s(2)=subplot(1,2,2);
        histogram(delta_through_bully,'BinWidth',.1)
        title('\Delta through tube chaser')
        xlabel('\Delta time through tube');
        ylabel('count');
        xlim([0 3]);
        set_fonts()
        linkaxes(s,'xy');
        
        f.Units = 'centimeters';
        f.Position = [3 3 8 4];
        exportgraphics(f,fullfile(ops.save_dir,['time_diff_through_tube',num2str(tt),'.png']));
    end
    
    % filter out following events above threshold
    close_following_index = close_following_index(delta_t_entry<ops.threshold_lag_at_detector & delta_t_exit<ops.threshold_lag_at_detector ...
        & delta_through_victim<ops.threshold_lag_through_tube & delta_through_bully<ops.threshold_lag_through_tube);

    % filter above threshold
%     close_following_index = close_following_index(delta_through_victim<ops.threshold_lag_through_tube & delta_through_bully<ops.threshold_lag_through_tube);


    % parse to chase-variables that indexes the original filtered data and
    % not the tube data
    chasing_index = cat(1,chasing_index,find(ismember(filtered_detections.Detection_start,tube_filtered.Detection_start(close_following_index))));
end

chasing_index = sort(chasing_index);

% threshold for multiple chasings through tube to collapse into one event?
%     putative_collapse = find(diff(filtered_detections.Detection_start(chasing_index))<threshold_lag_between_tubes);
%     for ii=1:numel(putative_collapse)
%         filtered_detections.Animal(chasing_index(putative_collapse(ii)))
%         filtered_detections.Animal(chasing_index(putative_collapse(ii)+1))
%     end

% write the chasing events to a day-by-day match matrix like for the
% tube test
unique_days = unique(filtered_detections.Day);
for day = 1:numel(unique_days)
   cur_chasings = chasing_index(ismember(filtered_detections.Day(chasing_index),unique_days{day}));

   % initialize matrix
   full_hierarchy(day).match_matrix_chasing = zeros(size(full_hierarchy(day).match_matrix));
   full_hierarchy(day).match_info_chasing = cell(size(full_hierarchy(day).match_info));
   for ii = 1:numel(cur_chasings)
       loser = filtered_detections.Animal{cur_chasings(ii)}; % == being chased
       loser_id = find(contains(full_hierarchy(1).ID,loser));
       winner = filtered_detections.Animal{cur_chasings(ii)+1}; % == chasing
       winner_id = find(contains(full_hierarchy(1).ID,winner));
       assert(~strcmp(winner,loser));

       % parse
       full_hierarchy(day).match_matrix_chasing(winner_id, loser_id) = full_hierarchy(day).match_matrix_chasing(winner_id, loser_id)+1;
       num = full_hierarchy(day).match_matrix_chasing(winner_id, loser_id);
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).winner_entry_time = filtered_detections.Time_sec(cur_chasings(ii)+1);
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).loser_entry_time = filtered_detections.Time_sec(cur_chasings(ii));
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).winner_exit_time = filtered_detections.Time_sec(cur_chasings(ii)+3);
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).loser_exit_time = filtered_detections.Time_sec(cur_chasings(ii)+2);
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).winner_ID = winner;
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).loser_ID = loser;
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).winner_entry_Detector = filtered_detections.Detector(cur_chasings(ii)+1);
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).winner_exit_Detector = filtered_detections.Detector(cur_chasings(ii)+3);
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).loser_entry_Detector = filtered_detections.Detector(cur_chasings(ii));
       full_hierarchy(day).match_info_chasing{winner_id,loser_id}(num).loser_exit_Detector = filtered_detections.Detector(cur_chasings(ii)+2);
   end

   % plot single chasing events for visual inspection
%     if 1==1
%         cur_data=full_hierarchy(day);
%         cur_data.save_dir = fullfile(save_dir,'plots','curation_v2',['group_',num2str(gr)],['day_',num2str(day)]);
%         if ~isfolder(cur_data.save_dir); mkdir(cur_data.save_dir); end
%         plot_chasing_events(cur_data);
%         close all;
%     end
end
end
