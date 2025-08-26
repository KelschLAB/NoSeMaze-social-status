function plot_tube_chasings_daycycle(data_dir,save_dir,cohortsTbl)
%%
% get daytime of the chasing events
collect_daytime = [];

% Loop through each cohort listed in the cohorts table
for i = 1:height(cohortsTbl)
    
    % Extract the current cohort's name from the table
    cohort = cohortsTbl.cohort{i};

    % Load chasing data
    cohort_data_dir = fullfile(data_dir,'processed',cohort,'tube','full_hierarchy_files');
    load(fullfile(cohort_data_dir,['full_hierarchy_' cohort '.mat']));

    for day = 1:numel(full_hierarchy)
        [locx,locy] = find(full_hierarchy(day).match_matrix_chasing);
        for dyad = 1:numel(locx)
            for events = 1:numel(full_hierarchy(day).match_info_chasing{locx(dyad),locy(dyad)})
                collect_daytime = cat(1,collect_daytime,full_hierarchy(day).match_info_chasing{locx(dyad),locy(dyad)}(events).winner_entry_time/3600);
            end
        end
    end
    
end

%%
% Plot a simple polar histogram with all events from all animals
f = figure;
h = polarhistogram(collect_daytime*(2*pi/24),24,'Normalization','probability',...
    'FaceColor',[227 30 36]./255,'FaceAlpha',1);
h.EdgeColor = 'k'; %[227 30 36]./255;
h.LineWidth = 1;
axis = gca;
axis.ThetaTickLabel = cellfun(@num2str,num2cell((0:2:22)'),'UniformOutput',0);
axis.FontSize = 6;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [3 3 4.5 3]);

%%
% save
[annot, srcInfo] = docDataSrc(gcf,fullfile(save_dir),mfilename('fullpath'),logical(1));
exportgraphics(gcf, fullfile(save_dir,'tube_chasings_daycycle.pdf'),'ContentType','vector','BackgroundColor','none');

% export source data
writetable(table(h.BinCounts',h.Values','Variablenames',{'BinCounts','BinProp'}),fullfile(save_dir,'SourceData_tube_chasings_daycycle.xlsx'));

close all
