function lick_params = compute_reward_prediction(lickport_data,lick_params,plotPath)

% The task contains an inert structure: 
% The probability of a reward is lower when the trial before was rewarded.
% See below:
CSplus_trials = find([lickport_data.reward]==1);

lick_params.RewProb_CSplus_minus1 = sum([lickport_data(CSplus_trials(2:end)-1).reward]==1)/sum([lickport_data(CSplus_trials(2:end)-1).reward]==0);
lick_params.RewProb_CSplus_minus2 = sum([lickport_data(CSplus_trials(3:end)-2).reward]==1)/sum([lickport_data(CSplus_trials(3:end)-2).reward]==0);
lick_params.RewProb_CSplus_minus3 = sum([lickport_data(CSplus_trials(4:end)-3).reward]==1)/sum([lickport_data(CSplus_trials(4:end)-3).reward]==0);
lick_params.RewProb_CSplus_minus4 = sum([lickport_data(CSplus_trials(5:end)-4).reward]==1)/sum([lickport_data(CSplus_trials(5:end)-4).reward]==0);
lick_params.RewProb_CSplus_minus5 = sum([lickport_data(CSplus_trials(6:end)-5).reward]==1)/sum([lickport_data(CSplus_trials(6:end)-5).reward]==0);
lick_params.RewProb_CSplus_minus6 = sum([lickport_data(CSplus_trials(7:end)-6).reward]==1)/sum([lickport_data(CSplus_trials(7:end)-6).reward]==0);

% separation into rewarded/non-rewarded before
x = 2:length(lickport_data);
trial_before_rewarded = logical([lickport_data(x-1).reward]);
baseline_licks_last_trial_rewarded = {lickport_data(x(trial_before_rewarded)).licks_bef_od};
baseline_licks_last_trial_non_rewarded = {lickport_data(x(~trial_before_rewarded)).licks_bef_od};
baseline_licks_last_trial_rewarded_omitfirst = cellfun(@(x) x(x>0.05),baseline_licks_last_trial_rewarded,'UniformOutput',0);
baseline_licks_last_trial_non_rewarded_omitfirst = cellfun(@(x) x(x>0.05),baseline_licks_last_trial_non_rewarded,'UniformOutput',0);

% figure
fig(2)=figure('Visible','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.3]);
subplot(2,2,1);
plot(movmean(cell2mat(cellfun(@numel, baseline_licks_last_trial_rewarded_omitfirst,'UniformOutput',0)).*2,[0,50]));
hold on; plot(movmean(cell2mat(cellfun(@numel, baseline_licks_last_trial_non_rewarded_omitfirst,'UniformOutput',0)).*2,[0,50]));
ax=gca; box off; ax.XLabel.String={'trials'};  ax.YLabel.String={'lick rate','(Hz, movmean over 50 trials)'}; ll=legend({'last: rew','last: no rew'},'Interpreter','none');
        
% statistics
nLicks_lastRew = cell2mat(cellfun(@numel, baseline_licks_last_trial_rewarded_omitfirst,'UniformOutput',0));
nLicks_lastNonRew = cell2mat(cellfun(@numel, baseline_licks_last_trial_non_rewarded_omitfirst,'UniformOutput',0));
[p, observeddifference, effectsize] = permutationTest(nLicks_lastNonRew, nLicks_lastRew, 10000);

lick_params.baseline_last_trial_rewarded_omitfirst = mean(cell2mat(cellfun(@numel, baseline_licks_last_trial_rewarded_omitfirst,'UniformOutput',0)));
lick_params.baseline_last_trial_non_rewarded_omitfirst = mean(cell2mat(cellfun(@numel, baseline_licks_last_trial_non_rewarded_omitfirst,'UniformOutput',0)));
lick_params.prediction_normalized_lick_diff = (lick_params.baseline_last_trial_non_rewarded_omitfirst-lick_params.baseline_last_trial_rewarded_omitfirst)./(lick_params.baseline_last_trial_rewarded_omitfirst+lick_params.baseline_last_trial_non_rewarded_omitfirst);
lick_params.prediction_lick_diff = (lick_params.baseline_last_trial_non_rewarded_omitfirst-lick_params.baseline_last_trial_rewarded_omitfirst);
lick_params.prediction_effectsize = effectsize;
lick_params.prediction_significance = p;

title({['last rew.: ' num2str(lick_params.baseline_last_trial_rewarded_omitfirst) ', last non-rew.: ' num2str(lick_params.baseline_last_trial_non_rewarded_omitfirst)],['effectsize = ' num2str(effectsize)]});

% save
[annot, srcInfo] = docDataSrc(fig(2),fullfile(plotPath),mfilename('fullpath'),logical(1));
exportgraphics(fig(2),fullfile(plotPath,['Baseline_Licks_Movmean_' lick_params.ID '.pdf']),'Resolution',300);
exportgraphics(fig(2),fullfile(plotPath,['Baseline_Licks_Movmean_' lick_params.ID '.png']),'Resolution',300);

% write csv
writetable(table(movmean(cell2mat(cellfun(@numel, baseline_licks_last_trial_rewarded_omitfirst,'UniformOutput',0)).*2,[0,50])','VariableNames',{'Licks_Last_Rew'}),fullfile(plotPath,['baseline_licks_movmean_lastRew_' lick_params.ID '.csv']));
writetable(table(movmean(cell2mat(cellfun(@numel, baseline_licks_last_trial_non_rewarded_omitfirst,'UniformOutput',0)).*2,[0,50])','VariableNames',{'Licks_Last_Non_Rew'}),fullfile(plotPath,['baseline_licks_movmean_lastNonRew_' lick_params.ID '.csv']));



end
