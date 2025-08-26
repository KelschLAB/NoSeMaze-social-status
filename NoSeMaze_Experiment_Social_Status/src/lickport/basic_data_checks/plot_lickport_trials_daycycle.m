function plot_lickport_trials_daycycle(data_dir,save_dir,cohortsTbl,exclusion_table)
%%
% get daytime of the competition events
collect_daytime = [];

% Loop through each cohort listed in the cohorts table
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

    % Loop over animals
    for an = 1:size(d,1)
        % time string
        time_strs = {d(an).events.time_StartTrial};
        % Compute total time in seconds
        time_in_seconds = cellfun(@(t) sum(sscanf(t, '%d:%d:%f') .* [3600; 60; 1]), time_strs);
        %
        collect_daytime = [collect_daytime;(time_in_seconds./3600)'];
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

axis = gca;
axis.ThetaZeroLocation = 'top';
axis.ThetaDir = 'clockwise';
axis.ThetaTick = 0:30:330;                                % every 2 h
axis.ThetaTickLabel = compose('%02d', 0:2:22);            % 0..22 by 2
axis.FontSize = 8;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [3 3 4.5 3]);

%%
% save
[annot, srcInfo] = docDataSrc(gcf,fullfile(save_dir),mfilename('fullpath'),logical(1));
exportgraphics(gcf, fullfile(save_dir,'lickport_trials_daycycle.pdf'),'ContentType','vector','BackgroundColor','none');

% export source data
writetable(table(h.BinCounts',h.Values','Variablenames',{'BinCounts','BinProp'}),fullfile(save_dir,'SourceData_lickport_trials_daycycle.xlsx'));

close all
