function plot_PSTH_single_animals(d, plotPath, cohort)

% Loop over each animal's data
for an = 1:numel(d)
    % Extract lickport data for this animal
    lickport_data = d(an).events;

    % Store the animal ID for later use
    lick_params.ID = lickport_data(1).ID;

    % Identify reversal indices
    reversal_index = [1, find(diff([lickport_data.phase]) == 1) + 1];
    lick_params.reversal_index = reversal_index;

    % Define trials of interest
    % 1. before reversal
    cs_plus_trials_before_reversal = [];
    cs_minus_trials_before_reversal = [];
    % Find trials for CS+ and CS- before reversal
    for rev = 2:numel(reversal_index)
        cs_plus_trials_before_reversal = cat(1, cs_plus_trials_before_reversal, reversal_index(rev-1) - 1 + find([lickport_data(reversal_index(rev-1):reversal_index(rev)-1).reward] == 1, 150, 'last')');
        cs_minus_trials_before_reversal = cat(1, cs_minus_trials_before_reversal, reversal_index(rev-1) - 1 + find([lickport_data(reversal_index(rev-1):reversal_index(rev)-1).reward] == 0, 150, 'last')');
    end
    lick_params.cs_plus_trials_before_reversal = cs_plus_trials_before_reversal;
    lick_params.cs_minus_trials_before_reversal = cs_minus_trials_before_reversal;
    % 2. after reversal
    cs_plus_trials_after_reversal = [];
    cs_minus_trials_after_reversal = [];
    % Find trials for CS+ and CS- after reversal
    for rev = 2:numel(reversal_index)
        if rev < numel(reversal_index)
            cs_plus_trials_after_reversal = cat(1, cs_plus_trials_after_reversal, reversal_index(rev)-1 + ...
                find([lickport_data(reversal_index(rev):reversal_index(rev+1)-1).reward] == 1, 150, 'first')');
            cs_minus_trials_after_reversal = cat(1, cs_minus_trials_after_reversal, reversal_index(rev)-1 + ...
                find([lickport_data(reversal_index(rev):reversal_index(rev+1)-1).reward] == 0, 150, 'first')');
        else
            cs_plus_trials_after_reversal = cat(1, cs_plus_trials_after_reversal, reversal_index(rev)-1 + ...
                find([lickport_data(reversal_index(rev):end).reward] == 1, 150, 'first')');
            cs_minus_trials_after_reversal = cat(1, cs_minus_trials_after_reversal, reversal_index(rev)-1 + ...
                find([lickport_data(reversal_index(rev):end).reward] == 0, 150, 'first')');
        end
    end

    % PSTH plot
    % Define trial ranges:
    phases_before = unique([lickport_data(cs_plus_trials_before_reversal).phase]);
    range_start = find([lickport_data(cs_plus_trials_before_reversal).phase] == min(phases_before));
    range_end = find([lickport_data(cs_plus_trials_before_reversal).phase] == max(phases_before));
    range_all = 1:length(cs_plus_trials_before_reversal);
    trial_ranges_plus_before = {range_all,range_start,range_end};

    phases_before = unique([lickport_data(cs_minus_trials_before_reversal).phase]);
    range_start = find([lickport_data(cs_minus_trials_before_reversal).phase] == min(phases_before));
    range_end = find([lickport_data(cs_minus_trials_before_reversal).phase] == max(phases_before));
    range_all = 1:length(cs_minus_trials_before_reversal);
    trial_ranges_minus_before = {range_all,range_start,range_end};

    range_name_before = {'all reversals: <=150 trials before','reversal #1: <=150 trials before','last reversal: <=150 trials before'};

    phases_after = unique([lickport_data(cs_plus_trials_after_reversal).phase]);
    range_start = find([lickport_data(cs_plus_trials_after_reversal).phase] == min(phases_after));
    range_end = find([lickport_data(cs_plus_trials_after_reversal).phase] == max(phases_after));
    range_all = 1:length(cs_plus_trials_after_reversal);
    trial_ranges_plus_after = {range_all,range_start,range_end};

    phases_after = unique([lickport_data(cs_minus_trials_after_reversal).phase]);
    range_start = find([lickport_data(cs_minus_trials_after_reversal).phase] == min(phases_after));
    range_end = find([lickport_data(cs_minus_trials_after_reversal).phase] == max(phases_after));
    range_all = 1:length(cs_minus_trials_after_reversal);
    trial_ranges_minus_after = {range_all,range_start,range_end};

    range_name_after = {'all reversals: <=150 trials after','reversal #1: <=150 trials after','last reversal: <=150 trials after'};

    % Initiate Figure
    fig(1)=figure('Visible','off');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.8, 0.5]);
    subpl_counter = 1;
    % Loop over trial ranges
    for range_idx = 1:length(trial_ranges_plus_before)
        % Subplot
        subplot(2,3,subpl_counter)
        % Input
        licks_cs_plus = [[lickport_data(cs_plus_trials_before_reversal(trial_ranges_plus_before{range_idx})).licks_bef_od],[lickport_data(cs_plus_trials_before_reversal(trial_ranges_plus_before{range_idx})).licks_aft_od]];
        licks_cs_minus = [[lickport_data(cs_minus_trials_before_reversal(trial_ranges_minus_before{range_idx})).licks_bef_od],[lickport_data(cs_minus_trials_before_reversal(trial_ranges_minus_before{range_idx})).licks_aft_od]];
        edges = 0:0.05:3.5;
        histcount_cs_plus = histcounts(licks_cs_plus,edges).*(20/length(cs_plus_trials_before_reversal(trial_ranges_plus_before{range_idx})));
        histcount_cs_minus = histcounts(licks_cs_minus,edges).*(20/length(cs_minus_trials_before_reversal(trial_ranges_minus_before{range_idx})));
        plot_PSTH(histcount_cs_plus,histcount_cs_minus,edges,range_name_before{range_idx})
        % write csv
        writetable(table([1:length(histcount_cs_plus)]',histcount_cs_plus',histcount_cs_minus','VariableNames',{'Edge','CSplus','CSminus'}),fullfile(plotPath,['SourceData_PSTH_Histogramm' num2str(subpl_counter) '_' lick_params.ID '_' cohort '.csv']));
        % counter update
        subpl_counter = subpl_counter+1;
    end
    for range_idx = 1:length(trial_ranges_plus_after)
        % Subplot
        subplot(2,3,subpl_counter)
        % Input
        licks_cs_plus = [[lickport_data(cs_plus_trials_after_reversal(trial_ranges_plus_after{range_idx})).licks_bef_od],[lickport_data(cs_plus_trials_after_reversal(trial_ranges_plus_after{range_idx})).licks_aft_od]];
        licks_cs_minus = [[lickport_data(cs_minus_trials_after_reversal(trial_ranges_minus_after{range_idx})).licks_bef_od],[lickport_data(cs_minus_trials_after_reversal(trial_ranges_minus_after{range_idx})).licks_aft_od]];
        edges = 0:0.05:3.5;
        histcount_cs_plus = histcounts(licks_cs_plus,edges).*(20/length(cs_plus_trials_after_reversal(trial_ranges_plus_after{range_idx})));
        histcount_cs_minus = histcounts(licks_cs_minus,edges).*(20/length(cs_minus_trials_after_reversal(trial_ranges_minus_after{range_idx})));
        plot_PSTH(histcount_cs_plus,histcount_cs_minus,edges,range_name_after{range_idx})
        % write csv
        writetable(table([1:length(histcount_cs_plus)]',histcount_cs_plus',histcount_cs_minus','VariableNames',{'Edge','CSplus','CSminus'}),fullfile(plotPath,['SourceData_PSTH_Histogramm' num2str(subpl_counter) '_' lick_params.ID '_' cohort '.csv']));
        % counter update
        subpl_counter = subpl_counter+1;
    end
    %
    sgtitle(lick_params.ID);
    % save
    [~, ~] = docDataSrc(fig(1),fullfile(plotPath),mfilename('fullpath'),logical(1));
    exportgraphics(fig(1),fullfile(plotPath,['PSTH_' lick_params.ID '_' cohort '.pdf']),'Resolution',300);
    exportgraphics(fig(1),fullfile(plotPath,['PSTH_' lick_params.ID '_' cohort '.png']),'Resolution',300);

end
end

function plot_PSTH(histcount_cs_plus,histcount_cs_minus,edges,range_name)
histogram('BinEdges', edges, 'BinCounts', histcount_cs_plus,'EdgeColor','none','FaceColor',[0,.5,1]);
hold on;
histogram('BinEdges', edges, 'BinCounts', histcount_cs_minus,'EdgeColor','none','FaceColor',[1,.5,.5]);
axis square;
box off;
ax = gca;
ax.XLim=[0,3.5];
ax.XTick=[0:.5:3.5];
ax.XLabel.String = 'time (s)';
ax.YLabel.String = 'lick rate (Hz)';
ax.YLim=[0,15];
ax.YTick=[0:5:15];
ll=legend({'CS+','CS-'},'FontSize',10);
ax.FontSize = 10;
tt=title([range_name]);
end