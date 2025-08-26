function [lickport_trials_per_day] = plot_lickport_trials_n_per_day(data_dir,save_dir,cohortsTbl,exclusion_table)
%%
% Assign animal counter
animal_counter = 1;

% Loop over cohorts
for i = 1:height(cohortsTbl)
    % Clearing
    clear exclusion_vector d

    % Current cohort
    cohort = cohortsTbl.cohort{i};
    disp(['Processing ' cohort]);
   
    % Define path
    lickport_data_dir = fullfile(data_dir,'interim',cohort,'lickport');

    % Load data
    load(fullfile(lickport_data_dir,['reorganizedData_' cohort '.mat']));
    
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

    % Loop over animals to count trials
    for an = 1:size(d,1)
        days = unique({d(an).events.date});
        for dd = 1:numel(days)
            lickport_trials_per_day(animal_counter,dd) = nnz(contains({d(an).events.date},days{dd}));
        end
        animal_counter=animal_counter+1;
    end 
end

%%
% plot mean and sem total competition events per day over time
figure;
line([0 22],[mean(lickport_trials_per_day,'all','omitnan') mean(lickport_trials_per_day,'all','omitnan')],'LineStyle','--','Color',[.5 .5 .5])
hold on
iosr.statistics.boxPlot(lickport_trials_per_day(:,1:21),...
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
ylabel('number of competition events');
title([num2str(height(cohortsTbl)) ' cohorts']);
xlim([0 22])
set(gcf, 'Units', 'centimeters');
set_fonts();
set(gcf, 'Position', [3 3 9 5]);

%%
% save
[annot, srcInfo] = docDataSrc(gcf,fullfile(save_dir),mfilename('fullpath'),logical(1));
exportgraphics(gcf, fullfile(save_dir,'lickport_trials_n_per_day.pdf'),'ContentType','vector','BackgroundColor','none');

% export source data
writetable(array2table(lickport_trials_per_day),fullfile(save_dir,'SourceData_lickport_trials_n_per_day.xlsx'));

close all