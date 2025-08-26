function plot_tube_competitions(hierarchy_data)
%% This function visualizes tube test from autonomouse.
% 15 seconds before and after the entry timepoint of the winner ist plotted
% and colorcoded for the loser and also all other animals 
%   red: winner entry gate colorcode
%   green: winner exit gate color == loser gate
%   blue: other tube
%
%
%% parse input
match_matrix = hierarchy_data.match_matrix;
match_info = hierarchy_data.match_info;
Data = hierarchy_data.Data;

%% loop over single events and create visualization
% for match = 1:sum(match_matrix,'all')
for winner = 1:size(match_matrix,1)
    for loser = 1:size(match_matrix,1)
        if isempty(match_info{winner,loser}); continue; end
        for ii=1:size(match_info{winner,loser},2)
           f=plot_single_event(Data,match_info{winner,loser}(1,ii));
           if isfield(hierarchy_data,'save_dir')
               saveas(f,fullfile(hierarchy_data.save_dir,[hierarchy_data.ID{winner},'vs',hierarchy_data.ID{loser},'_',num2str(ii),'.png']));
           end
        end
    end
end
end


%%
function f=plot_single_event(Data,this_match_info)

time_base = round(this_match_info.winner_entry_time-15:.5:this_match_info.winner_exit_time+15,1);
% time_base = this_match_info.winner_entry_time-15:this_match_info.winner_exit_time+15;
relevant_Data = Data(Data.Time_sec>time_base(1) & Data.Time_sec<time_base(end),:);

winner_heat = zeros(1,numel(time_base));
% winner_heat(any(time_base==round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, this_match_info.winner_ID) & relevant_Data.Detector == this_match_info.winner_entry_Detector),1),1)) = 1; % entry detector
% winner_heat(any(time_base==round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, this_match_info.winner_ID) & relevant_Data.Detector == this_match_info.winner_exit_Detector),1),1)) = 2; % exit detector

% sub-second temporal resolution
winner_heat(any(histcounts(round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, this_match_info.winner_ID) & relevant_Data.Detector == this_match_info.winner_entry_Detector),1),time_base),1)) = 1; % entry detector
winner_heat(any(histcounts(round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, this_match_info.winner_ID) & relevant_Data.Detector == this_match_info.winner_exit_Detector),1),time_base),1)) = 2; % exit detector

loser_heat = zeros(1,numel(time_base));
% loser_heat(any(time_base==round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, this_match_info.loser_ID) & relevant_Data.Detector == this_match_info.winner_entry_Detector),1),1)) = 1; % entry detector
% loser_heat(any(time_base==round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, this_match_info.loser_ID) & relevant_Data.Detector == this_match_info.winner_exit_Detector),1),1)) = 2; % exit detector
loser_heat(any(histcounts(round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, this_match_info.loser_ID) & relevant_Data.Detector == this_match_info.winner_entry_Detector),1),time_base),1)) = 1; % entry detector
loser_heat(any(histcounts(round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, this_match_info.loser_ID) & relevant_Data.Detector == this_match_info.winner_exit_Detector),1),time_base),1)) = 2; % exit detector




% other animals
unique_ids = unique(Data.Animal);
unique_ids(ismember(unique_ids,{this_match_info.winner_ID,this_match_info.loser_ID}))=[];
others_heat = zeros(numel(unique_ids),numel(time_base));
for an = 1:size(others_heat,1)
%     others_heat(an, any(time_base==round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, unique_ids{an}) & relevant_Data.Detector == this_match_info.winner_entry_Detector)))) = 1; % entry detector
%     others_heat(an, any(time_base==round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, unique_ids{an}) & relevant_Data.Detector == this_match_info.winner_exit_Detector)))) = 2; % exit detector
%     others_heat(an, any(time_base==round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, unique_ids{an}) & relevant_Data.Detector ~= this_match_info.winner_exit_Detector ...
%         & relevant_Data.Detector ~= this_match_info.winner_entry_Detector)))) = 3; % other detector

    % sub-second temporal resolution
    others_heat(an, any(histcounts(round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, unique_ids{an}) & relevant_Data.Detector == this_match_info.winner_entry_Detector),1),time_base),1)) = 1; % entry detector
    others_heat(an, any(histcounts(round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, unique_ids{an}) & relevant_Data.Detector == this_match_info.winner_exit_Detector),1),time_base),1)) = 2; % exit detector
    others_heat(an, any(histcounts(round(relevant_Data.Time_sec(ismember(relevant_Data.Animal, unique_ids{an}) & relevant_Data.Detector ~= this_match_info.winner_exit_Detector ...
        & relevant_Data.Detector ~= this_match_info.winner_entry_Detector),1),time_base),1)) = 3; % other detector
end

f=figure('Visible','off');
subplot(10,1,[1 3]);
imagesc(time_base,1:2,cat(1,winner_heat,loser_heat),[0 3])
yticks([1,2]);
yticklabels({this_match_info.winner_ID,this_match_info.loser_ID})
xlabel('seconds');
title({['event ID: ', num2str(this_match_info.unique_event_ID),'. multiple-animals: ', num2str(this_match_info.multiple_animals)],['entry: ' num2str(this_match_info.winner_entry_Detector) ', exit: ' num2str(this_match_info.winner_exit_Detector)]});
set_fonts()

subplot(10,1,[5 10]);
imagesc(time_base,1:size(others_heat,1),others_heat,[0 3])
yticks(1:size(others_heat,1));
yticklabels(unique_ids)
xlabel('seconds');
set_fonts()

map = [1 1 1
    1 0 0
    0 1 0
    0 0 1];
colormap(map);

end

%%
function set_fonts(axis)

if nargin < 1
    axis = gca;
end

set(get(axis, 'XLabel'), 'FontSize', 9);%6);
set(get(axis, 'YLabel'), 'FontSize', 9);%6);
set(axis, 'FontSize', 9);%6);
set(axis, 'FontName', 'Arial');

set(get(axis, 'Title'), 'FontSize', 8);%8);


end