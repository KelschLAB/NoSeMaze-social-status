%% analyse_association_group_level.m
% Author: Jonathan Reinwald
% Date: July 2025




% Pre-clearing
clear all; clc; close all;

% Set pathes
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory
output_dir = fullfile(main_dir,'results','figures','cross_cohort','associations_group_level');
processed_dir = fullfile(main_dir, 'data', 'processed');
if ~isfolder(output_dir)
    mkdir(output_dir);
end
addpath(genpath(fullfile(main_dir,'src','helpers')));

%% Load Summary Data
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');
if ~isfile(summary_file)
    error('Summary file not found. Please run `main_build_cross_cohort_summary` first.');
end
load(summary_file, 'summary_data');

%% Load cohort characteristics (cohort-level information on e.g., transitivity and so on)
CohortCharacteristics = readtable(fullfile(main_dir,'data','processed','cross_cohort_files','tube','CohortCharacteristics_TubeCompetitions_D1_21_Koptimal.csv'));

%% Cohort Selection
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));
% select cohorts
cohortsTbl = cohortsTbl((cohortsTbl.use_tube == 1), :);

%% Filter input data (cohort selection)
fields = fieldnames(summary_data.tube);
for field_idx = 1:numel(fields)
    summary_data.tube.(fields{field_idx}) = summary_data.tube.(fields{field_idx})(ismember(summary_data.tube.(fields{field_idx}).cohort,cohortsTbl.cohort),:);
end

%% Define relevant variables
DSz_Competition = summary_data.tube.D1_21.DSz_Competition;
Fraction_Wins_Chasing = summary_data.tube.D1_21.Fraction_Wins_Chasing;
Fraction_Losses_Chasing = summary_data.tube.D1_21.Fraction_Losses_Chasing;
N_Events_Chasing = summary_data.tube.D1_21.N_Events_Chasing;
cohorts = summary_data.tube.D1_21.cohort;
% Identify unique cohorts
myCohorts = unique(cohorts);

%% Compute within-cohort correlations
corrCoeff_ChasingToCompetition = nan(length(myCohorts),1);
corrCoeff_ChasingToCompetition_pearson = nan(length(myCohorts),1);
sum_ChasingEvents = nan(length(myCohorts),1);
FrChas_TopRankedAnimal = nan(length(myCohorts),1);

% figure
fig = figure('Visible', 'on');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.9]);

for ix = 1:length(myCohorts)
    mask = strcmp(cohorts, myCohorts{ix});
    
    % Spearman correlation: Fraction_Wins_Chasing vs. DSz_Competition
    corrCoeff_ChasingToCompetition(ix,1) = corr(Fraction_Wins_Chasing(mask), DSz_Competition(mask), 'type', 'Spearman');
    corrCoeff_ChasingToCompetition_pearson(ix,1) = corr(Fraction_Wins_Chasing(mask), DSz_Competition(mask), 'type', 'Pearson');

    % N_Events per cohort
    sum_ChasingEvents(ix,1) = sum(N_Events_Chasing(mask));

    [a,b]=sort(DSz_Competition(mask),'descend');
    FractionChasingWinsTemp = Fraction_Wins_Chasing(mask);
    FrChas_TopRankedAnimal(ix,1) = nansum(FractionChasingWinsTemp(b(1)));
    FractionChasingLossesTemp = Fraction_Losses_Chasing(mask);
    FrChasLosses_TopRankedAnimal(ix,1) = nansum(FractionChasingLossesTemp(b(1)));
    ChasingEventsTemp = N_Events_Chasing(mask);
    NChas_TopRankedAnimal(ix,1) = sum(ChasingEventsTemp(b(1)))./sum(ChasingEventsTemp);
    %
    subplot(4,5,ix)
    scatter(DSz_Competition(mask),Fraction_Wins_Chasing(mask),'filled','MarkerFaceColor',[0.8,0.2,0.2],'MarkerEdgeColor','none','Marker','o','SizeData',20);
    hold on;
    ylabel('Fraction Wins (chasing)','Interpreter','none');
    xlabel('DSz (tube)','Interpreter','none');
    xlim([-2.2,2.2]); 
    ll=lsline; ll.Color=[0.8,0.5,0.5];
    ylim([0,max(Fraction_Wins_Chasing(mask))*1.1]); yy = ylim;% max(Fraction_Wins_Chasing(mask))
    tx=text(-1.5,0.9*yy(2),['r_S_p_e_a_r_m_a_n = ' num2str(round(corrCoeff_ChasingToCompetition(ix,1),2))]);
    tx=text(-1.5,0.8*yy(2),['r_P_e_a_r_s_o_n = ' num2str(round(corrCoeff_ChasingToCompetition_pearson(ix,1),2))]);
    tx=text(-1.5,0.7*yy(2),['transitivity = ' num2str(round(CohortCharacteristics.transitivity_pt(ix),2))]);
    title(myCohorts{ix});
end

[annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(output_dir, sprintf('Scatter_Plots_Chasing_To_Competition_by_Groups.pdf')));

%% Merge correlation results into CohortCharacteristics table
% Match cohort names to ensure correct alignment
[isCohort, idxMap] = ismember(CohortCharacteristics.cohort, myCohorts);

% Add correlation coefficients to the table
CohortCharacteristics.corrCoeff_ChasingToCompetition = nan(height(CohortCharacteristics),1);
CohortCharacteristics.sum_ChasingEvents = nan(height(CohortCharacteristics),1);

CohortCharacteristics.corrCoeff_ChasingToCompetition(isCohort) = corrCoeff_ChasingToCompetition(idxMap(isCohort));
CohortCharacteristics.sum_ChasingEvents(isCohort) = sum_ChasingEvents(idxMap(isCohort));
CohortCharacteristics.FrChas_TopRankedAnimal(isCohort) = FrChas_TopRankedAnimal(idxMap(isCohort));

%% Load and integrate the quadrant quantification
% from: myRootPath\NoSeMaze_Experiment_Socia_Status\analysis\scripts\association_between_metrics\analyse_match_matrices.m
load('myRootPath\NoSeMaze_Experiment_Socia_Status\results\figures\cross_cohort\associations_between_metrics\tube\match_matrices\Quadrant_Quantification_ChasingSortedByComp_D1_21.mat')
load('myRootPath\NoSeMaze_Experiment_Socia_Status\results\figures\cross_cohort\associations_between_metrics\tube\match_matrices\ChasingDirection_D1_21.mat')
CohortCharacteristics.ChasingbyComp_TopLeftQuadrant = quadrant_quantification.chasing_by_comp_TopLeft;
CohortCharacteristics.ChasingbyComp_TopRightQuadrant = quadrant_quantification.chasing_by_comp_TopRight;
CohortCharacteristics.ChasingbyComp_BottomLeftQuadrant = quadrant_quantification.chasing_by_comp_BottomLeft;
CohortCharacteristics.ChasingbyComp_BottomRightQuadrant = quadrant_quantification.chasing_by_comp_BottomRight;
CohortCharacteristics.down_vs_up_chasing = chasing_direction.down_vs_up_chasing;
CohortCharacteristics.down_vs_up_chasing_weighted = chasing_direction.down_vs_up_chasing_weighted;
CohortCharacteristics.top_down_chasing = chasing_direction.top_down_chasing;
CohortCharacteristics.top_down_chasing_weighted = chasing_direction.top_down_chasing_weighted;
CohortCharacteristics.bottom_up_chasing = chasing_direction.bottom_up_chasing;
CohortCharacteristics.bottom_up_chasing_weighted = chasing_direction.bottom_up_chasing_weighted;

%% Correlation matrix of cohort-level characteristics
% Select variables of interest
vars = {'steepness', 'transitivity_tri', 'transitivity_pt', ...
        'linearity_h1', 'linearity_h2', 'stabilityIndex', ...
        'corrCoeff_ChasingToCompetition','sum_ChasingEvents', 'FrChas_TopRankedAnimal', ...
        'ChasingbyComp_TopLeftQuadrant','ChasingbyComp_TopRightQuadrant','ChasingbyComp_BottomLeftQuadrant','ChasingbyComp_BottomRightQuadrant', ...
        'down_vs_up_chasing','down_vs_up_chasing_weighted','top_down_chasing','top_down_chasing_weighted','bottom_up_chasing','bottom_up_chasing_weighted'};       

dataMatrix = table2array(CohortCharacteristics(:, vars));
[rmat, pmat] = corr(dataMatrix, 'type', 'Spearman');

%% Visualization
% figure
fig = figure('Visible', 'on');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.9]);

subplot(2,2,1);
imagesc(rmat);
axis square;
xticks(1:length(vars)); xticklabels(vars); xtickangle(45);
yticks(1:length(vars)); yticklabels(vars);
ax = gca;
ax.TickLabelInterpreter = 'none';
colorbar;
title('Spearman Correlation Matrix of Cohort Characteristics');
% Overlay significance markers
hold on;
for i = 1:size(rmat,1)
    for j = 1:size(rmat,2)
        if pmat(i, j) < (0.05)%pBonferroni)
            % Display both r and asterisk for significant values
            text(j, i, sprintf('%.2f*', rmat(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'white', 'FontSize', 8, 'FontWeight', 'bold');
        else
            % Display r without asterisk for non-significant values
            text(j, i, sprintf('%.2f', rmat(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'black', 'FontSize', 4, 'FontWeight', 'normal');
        end
    end
end

subplot(2,2,2);

% Extract data
x = CohortCharacteristics.corrCoeff_ChasingToCompetition;
y = CohortCharacteristics.transitivity_pt;

% Scatter plot
scatter(x, y, 'filled', ...
    'MarkerFaceColor', [0.8, 0.2, 0.2], ...
    'MarkerEdgeColor', 'none', ...
    'Marker', 'd', ...
    'SizeData', 20);
hold on;
axis square;
xlabel({'association strength (Spearmans r)', 'chasing to tube rank'}, 'Interpreter', 'none');
ylabel('transitivity', 'Interpreter', 'none');
xlim([-0.2, 0.8]);
ylim([0.7, 1]);
xticklabels([-0.2:0.2:0.8]);
yy = ylim;

% Add least-squares line
ll = lsline;
ll.Color = [0.8, 0.5, 0.5];

% Correlation text
tx = text(0.75, yy(2) - 0.1 * diff(yy), ...
    ['r_S_p_e_a_r_m_a_n = ' num2str(round(corr(y, x, 'type', 'Spearman'), 2))]);
tx = text(0.75, yy(2) - 0.2 * diff(yy), ...
    ['r_P_e_a_r_s_o_n = ' num2str(round(corr(y, x, 'type', 'Pearson'), 2))]);

% Save source data
writetable(array2table([x,y],'VariableNames',{'corrCoeff_ChasingToCompetition','transitivity'}), fullfile(output_dir, sprintf('SourceData_Transitivity_CorrCoeffChasComp.csv')));

%% Rank plot
rx = tiedrank(x);
ry = tiedrank(y);

subplot(2,2,3);
scatter(rx, ry, 'filled', ...
    'MarkerFaceColor', [0.8, 0.2, 0.2], ...
    'MarkerEdgeColor', 'none', 'SizeData', 20);
axis square;
xlabel('Rank of association strength');
ylabel('Rank of transitivity');
xlim([0, max(rx)+1]); ylim([0, max(ry)+1]);

% Optional: add regression line for visualization only
p = polyfit(rx, ry, 1);
hold on;
plot(rx, polyval(p, rx), 'Color', [0.8, 0.5, 0.5], 'LineWidth', 1.5);

% Add Spearman rho
yy = ylim;
text(max(rx)-4, yy(2)-0.1*diff(yy), ...
    ['r_{Spearman} = ' num2str(round(corr(x, y, 'type', 'Spearman'), 2))]);

% Fit linear model for CI
mdl = fitlm(rx, ry);

% Generate predictions and CI band
xpred = linspace(min(xlim), max(xlim), 100)';
[ypred, yCI] = predict(mdl, xpred);

% Plot CI band
fill([xpred; flipud(xpred)], [yCI(:,1); flipud(yCI(:,2))], ...
    [0.8, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

subplot(2,2,4);

% Extract data
x = CohortCharacteristics.FrChas_TopRankedAnimal;
y = CohortCharacteristics.transitivity_pt;

% Scatter plot
scatter(x, y, 'filled', ...
    'MarkerFaceColor', [0.8, 0.2, 0.2], ...
    'MarkerEdgeColor', 'none', ...
    'Marker', 'd', ...
    'SizeData', 20);
hold on;
axis square;
xlabel({'fraction of chasing of top-rank animal'}, 'Interpreter', 'none');
ylabel('transitivity', 'Interpreter', 'none');
% Add least-squares line
ll = lsline;
ll.Color = [0.8, 0.5, 0.5];

xlim([0,0.4]);
ylim([0.7, 1]);
xticklabels([0:0.1:0.4]);
yy = ylim;

% Correlation text
tx = text(0.05, yy(2) - 0.1 * diff(yy), ...
    ['r_S_p_e_a_r_m_a_n = ' num2str(round(corr(y, x, 'type', 'Spearman'), 2))]);
tx = text(0.05, yy(2) - 0.2 * diff(yy), ...
    ['r_P_e_a_r_s_o_n = ' num2str(round(corr(y, x, 'type', 'Pearson'), 2))]);

% Save source data
writetable(array2table([x,y],'VariableNames',{'FrChas_TopRankedAnimal','transitivity'}), fullfile(output_dir, sprintf('SourceData_Transitivity_FrChasTopRankAnimal.csv')));

% Plot figure
[annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(output_dir, sprintf('Association_Matrix_Cohort_Characteristics_to_Social_Behavior.pdf')));

%% Optional: Save updated table
writetable(CohortCharacteristics, fullfile(output_dir, 'CohortCharacteristics_WithCorrelations.csv'));
