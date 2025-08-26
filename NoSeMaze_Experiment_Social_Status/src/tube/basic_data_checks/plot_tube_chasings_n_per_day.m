function [chasings_per_day]=plot_tube_chasings_n_per_day(data_dir,save_dir,cohortsTbl)
%%
% Preallocate empty matrix
chasings_per_day = nan(height(cohortsTbl),30); % preallocate matrix to store cohorts x day number of matches

% Counter
cohort_counter=1;

% Loop through each cohort listed in the cohorts table
for i = 1:height(cohortsTbl)
    
    % Extract the current cohort's name from the table
    cohort = cohortsTbl.cohort{i};

    % Load chasing data
    cohort_data_dir = fullfile(data_dir,'processed',cohort,'tube','full_hierarchy_files');
    load(fullfile(cohort_data_dir,['full_hierarchy_' cohort '.mat']));

    % number of matches over days global
    total_chasings_per_day = [];
    for day = 1:numel(full_hierarchy)
        total_chasings_per_day = cat(1, total_chasings_per_day, sum(full_hierarchy(day).match_matrix_chasing,'all'));        
    end
    chasings_per_day(cohort_counter,1:numel(full_hierarchy)) = total_chasings_per_day;
    
    % update counter
    cohort_counter=cohort_counter+1;  
end

%%
% plot mean and sem total chasing events per day over time
figure;
line([0 22],[mean(chasings_per_day,'all','omitnan') mean(chasings_per_day,'all','omitnan')],'LineStyle','--','Color',[.5 .5 .5])
hold on
iosr.statistics.boxPlot(chasings_per_day(:,1:21),...
    'theme','colorboxes',...
    'boxColor',[0 0 0],...
    'showMean',false,...
    'meanMarker','d',...
    'meanColor',[1 1 1], ...
    'meanSize',4,...
    'medianColor',[1 1 1], ...
    'outlierSize',4,...
    'symbolMarker','o');
xlabel('days');
ylabel('number of chasing events');
title([num2str(height(cohortsTbl)) ' cohorts']);
xlim([0 22])
set(gcf, 'Units', 'centimeters');
set_fonts();
set(gcf, 'Position', [3 3 9 5]);

%%
% save
[annot, srcInfo] = docDataSrc(gcf,fullfile(save_dir),mfilename('fullpath'),logical(1));
exportgraphics(gcf, fullfile(save_dir,'tube_chasings_n_per_day.pdf'),'ContentType','vector','BackgroundColor','none');

% export source data
writetable(array2table(chasings_per_day),fullfile(save_dir,'SourceData_tube_chasings_n_per_day.xlsx'));

close all