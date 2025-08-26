%% Pre-Clearing
clear; clc; close all;

%% Set Directories
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory
processed_dir = fullfile(main_dir, 'data', 'processed');
results_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort', 'associations_between_metrics');
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');
plot_dir = fullfile(results_dir,'tube','LME_Chasing_Vs_Competition');
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

% Add required paths
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));
addpath(genpath(fullfile(main_dir, 'src', 'analysis', 'associations_between_metrics')));

%% Load Summary Data
if ~isfile(summary_file)
    error('Summary file not found. Please run `main_build_cross_cohort_summary` first.');
end
load(summary_file, 'summary_data');

%% Cohort Selection
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));
% Filtering
cohortsTbl = cohortsTbl(cohortsTbl.use_tube == 1, :);

% Filter input data (cohort selection)
fields = fieldnames(summary_data.tube);
for field_idx = 1:numel(fields)
    summary_data.tube.(fields{field_idx}) = summary_data.tube.(fields{field_idx})(ismember(summary_data.tube.(fields{field_idx}).cohort,cohortsTbl.cohort),:);
end

%% Create concatenated data table (competition and chasing)
DSz = [summary_data.tube.D1_End.DSz_Competition;summary_data.tube.D1_End.DSz_Chasing];
Wins = [summary_data.tube.D1_End.Fraction_Wins_Competition;summary_data.tube.D1_End.Fraction_Wins_Chasing];
Losses = [summary_data.tube.D1_End.Fraction_Losses_Competition;summary_data.tube.D1_End.Fraction_Losses_Chasing];
% Wins = boxcox(Wins+0.0001);
% Losses = boxcox(Losses+0.0001);
N_Events = [summary_data.tube.D1_End.N_Events_Competition;summary_data.tube.D1_End.N_Events_Chasing];
N_Events = boxcox(N_Events);
n = height(summary_data.tube.D1_End.Fraction_Wins_Competition);
Type = [repmat({'Competition'}, n, 1); repmat({'Chasing'}, n, 1)];
ID = [summary_data.tube.D1_End.Mouse_RFID;summary_data.tube.D1_End.Mouse_RFID];
myTable = table(DSz,Wins,Losses,N_Events,Type,ID,'VariableNames',{'DSz','Wins','Losses','N_Events','Type','ID'});
myTable = rmmissing(myTable);


% Fit the model allowing for random slopes on N_Events
lme_random_slope_and_intercept = fitlme(myTable, 'Wins ~ Losses * Type + (Losses * Type|ID)');
% Fit the mixed-effects model again
lme_random_intercept = fitlme(myTable, 'Wins ~ Losses * Type + (1|ID)');

% Prepare figure
fig = figure('Name', 'Comparison: Wins to Losses', ...
             'Color', 'w', ...
             'Visible', 'on', ...
             'Position', [100, 100, 1200, 800]);  % [left, bottom, width, height]

% Plot histogram of residuals (normality check)
subplot(2,3,1);
histogram(residuals(lme_random_slope_and_intercept), 30, 'Normalization', 'pdf');
hold on;
x = linspace(min(residuals(lme_random_slope_and_intercept)), max(residuals(lme_random_slope_and_intercept)), 100);
plot(x, normpdf(x, mean(residuals(lme_random_slope_and_intercept)), std(residuals(lme_random_slope_and_intercept))), 'r', 'LineWidth', 2);
xlabel('Residuals'); ylabel('Density');
title('Histogram of Residuals with Normal Fit');
hold off;

% Q-Q Plot
subplot(2,3,2);
qqplot(residuals(lme_random_slope_and_intercept));
title('Q-Q Plot of Residuals');

% Extract residuals
resid_vals = residuals(lme_random_slope_and_intercept);

% Lilliefors test (Kolmogorov-Smirnov for large samples)
[h_lillie, p_lillie] = lillietest(resid_vals);
disp(['Lilliefors Test p-value: ', num2str(p_lillie)]);

% Anderson-Darling test
[h_ad, p_ad] = adtest(resid_vals);
disp(['Anderson-Darling Test p-value: ', num2str(p_ad)]);

% Shapiro-Wilk test (not built-in, requires Statistics Toolbox)
try
    [h_sw, p_sw] = swtest(resid_vals);
    disp(['Shapiro-Wilk Test p-value: ', num2str(p_sw)]);
catch
    disp('Shapiro-Wilk Test not available in this MATLAB version.');
end

% Plot residuals vs. fitted values (homoscedasticity check)
subplot(2,3,3);
scatter(fitted(lme_random_intercept), residuals(lme_random_intercept));
xlabel('Fitted Values'); ylabel('Residuals');
title('Residuals vs. Fitted Values');
yline(0, '--r');

% Compare the models
disp('Random Slope + Intercept Model:'); disp(lme_random_slope_and_intercept);
disp('Random Intercept Model:'); disp(lme_random_intercept);

% Scatter plot with regression lines for Chasing and Competition
subplot(2,3,4);
gscatter(myTable.Losses, myTable.Wins, myTable.Type, 'br', 'os', 'filled');
hold on;
x = linspace(min(myTable.Losses), max(myTable.Losses), 100);
% Get fixed effects
FE = fixedEffects(lme_random_intercept);
x = linspace(min(myTable.Losses), max(myTable.Losses), 100);

% Competition (Type=0): intercept + slope * x
y_comp = FE(1) + FE(2) * x;

% Chasing (Type=1): intercept + type_effect + (losses + interaction)*x
y_chase = FE(1) + FE(3) + (FE(2) + FE(4)) * x;

% Plot
plot(x, y_comp, 'b-', 'LineWidth', 2);
plot(x, y_chase, 'r-', 'LineWidth', 2);
legend('Competition', 'Chasing', 'Competition Fit', 'Chasing Fit');
title('Wins vs. Losses Interaction');
xlabel('Losses'); ylabel('Wins');
grid on;
hold off;

% Source data
writetable(myTable,fullfile(plot_dir, sprintf('Comparison_LME_Chasing_Vs_Competition_WinsToLosses.csv')));

% Display AIC and BIC for both models
disp('Random Slope + Intercept Model (AIC, BIC):');
disp(lme_random_slope_and_intercept.ModelCriterion.AIC);
disp(lme_random_slope_and_intercept.ModelCriterion.BIC);
disp('Random Intercept Model (AIC, BIC):');
disp(lme_random_intercept.ModelCriterion.AIC);
disp(lme_random_intercept.ModelCriterion.BIC);

stat = compare(lme_random_intercept,lme_random_slope_and_intercept);
% Display the p-value
disp('Likelihood Ratio Test p-value:');
disp(stat.pValue);

disp('Likelihood Ratio Test result:');
disp(stat);

%% === Bootstrapping Fixed Effects: Cluster-Resampled Confidence Intervals ===

% Justification:
% Mixed-effects models assume normality and homoscedasticity, which were violated (see residual diagnostics).
% Therefore, we used a non-parametric bootstrap approach (resampling by subject) to assess the robustness of fixed effects.

nBoot = 100;
coef_names = {'Intercept','Losses','Type_Chasing','Losses:Type_Chasing'};
formula = 'Wins ~ Losses*Type + (1|ID)';
% formula = 'Wins ~ Losses*Type +(Losses * Type|ID)';
uniqueIDs = unique(myTable.ID);
nIDs = numel(uniqueIDs);
nCoef = numel(coef_names);
bootstat = nan(nBoot, nCoef);  % store fixed effects from each bootstrap

rng(42);  % for reproducibility

for b = 1:nBoot
    disp(num2str(b));
    resampledIDs = datasample(uniqueIDs, nIDs, 'Replace', true);
    bootData = [];

    for j = 1:length(resampledIDs)
        thisID = resampledIDs(j);
        if isstring(thisID)
            theseRows = myTable(myTable.ID == thisID, :);
        else
            theseRows = myTable(strcmp(myTable.ID, thisID), :);
        end
        bootData = [bootData; theseRows];  % Preserve duplication
    end

    try
        tempModel = fitlme(bootData, formula);
        bootstat(b, :) = fixedEffects(tempModel)';
    catch
        continue;
    end
end

% Remove failed fits
bootstat = bootstat(~any(isnan(bootstat),2), :);

% Compute bootstrap confidence intervals and SEs
ci_lower = quantile(bootstat, 0.025);
ci_upper = quantile(bootstat, 0.975);
boot_means = mean(bootstat);
boot_se = std(bootstat);

% Display
fprintf('\nBOOTSTRAP 95%% CONFIDENCE INTERVALS:\n');
for i = 1:numel(coef_names)
    fprintf('%s: %.4f ± %.4f (CI: [%.4f, %.4f])\n', ...
        coef_names{i}, boot_means(i), boot_se(i), ci_lower(i), ci_upper(i));
end

%% === Optional: Visualize Bootstrap CIs ===
subplot(2,3,6);
errorbar(1:nCoef, boot_means, boot_means - ci_lower, ci_upper - boot_means, ...
    'o', 'LineWidth', 2, 'CapSize', 8);
xlim([0.5 nCoef+0.5]);
xticks(1:nCoef);
xticklabels(coef_names);
ylabel('Fixed Effect Estimate');
xlabel('Model Term');
title('Bootstrapped 95% Confidence Intervals');
yline(0, '--k');
grid on;


writetable(array2table([boot_means;ci_lower;ci_upper],'Variablenames',{'intercept','losses','chasing','interaction'}),fullfile(plot_dir, sprintf('ErrorBarPlot.csv')));


% --- Save figure ---
[~, ~] = docDataSrc(fig,fullfile(plot_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(plot_dir, 'Comparison_LME_Chasing_Vs_Competition_WinsToLosses.pdf'));
disp(['Saved figure to: ' fullfile(plot_dir, 'Comparison_LME_Chasing_Vs_Competition_WinsToLosses.pdf')]);
