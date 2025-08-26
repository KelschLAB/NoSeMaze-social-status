function check_LME_assumptions(lme_model, varargin)
% check_LME_assumptions - Check residual diagnostics for a fitted LME model
%
% Usage:
%   check_LME_assumptions(lme_model)
%   check_LME_assumptions(lme_model, 'plot_save_path', 'path/to/save', 'fig_name', 'MyModel', 'run_robust', true)
%
% Inputs:
%   lme_model - fitted linear mixed-effects model (fitlme output)
% Optional Name-Value:
%   'plot_save_path' - full file path (without extension) to save the plot
%   'fig_name' - figure title name (string)
%   'run_robust' - true/false, whether to run a robust alternative model (default: false)

%% Parse Inputs
p = inputParser;
addOptional(p, 'plot_save_path', '');
addOptional(p, 'fig_name', 'LME Model Diagnostics');
addOptional(p, 'run_robust', false);
parse(p, varargin{:});
plot_save_path = p.Results.plot_save_path;
fig_name = p.Results.fig_name;
run_robust = p.Results.run_robust;

%% Step 1: Prepare residuals and fitted values
res = residuals(lme_model);
fit_vals = fitted(lme_model);

%% Step 2: Quantitative Residual Slope Test
mdl_res = fitlm(fit_vals, res);
pval_slope = mdl_res.Coefficients.pValue(2);
slope_estimate = mdl_res.Coefficients.Estimate(2);

%% Step 3: Plot residual diagnostics
fig = figure('Name', fig_name, 'Position', [100 100 1200 800], 'Visible','off');

subplot(2,2,1);
plot(fit_vals, res, 'ko', 'MarkerFaceColor', [0.6 0.6 0.6]);
xlabel('Fitted Values');
ylabel('Residuals');
title('Residuals vs Fitted Values');
grid on; box off; axis square;
lsline; % Least squares fit line

subplot(2,2,2);
histogram(res, 'Normalization', 'pdf');
xlabel('Residual');
ylabel('Density');
title('Histogram of Residuals');
grid on; box off; axis square;

subplot(2,2,3);
qqplot(res);
title('QQ Plot of Residuals');
grid on; box off; axis square;

subplot(2,2,4);
scatter(fit_vals, lme_model.Variables.Response, 'ko', 'filled');
xlabel('Fitted Values');
ylabel('Observed Response');
title('Observed vs Fitted');
grid on; box off; axis square;
lsline;

sg = sgtitle(fig_name,'FontWeight','bold');
sg.Interpreter = 'none';

%% Step 4: Create summary text
summary_text = {
    sprintf('Slope p = %.4g', pval_slope)
    sprintf('Slope = %.3f', slope_estimate)
    };
if pval_slope < 0.05
    summary_text = [summary_text; {'\bf⚠️ Trend detected!'}];
    box_color = [1 0.8 0.8]; % Light red
else
    summary_text = [summary_text; {'\bf✅ No trend'}];
    box_color = [0.8 1 0.8]; % Light green
end

%% Step 5: (Optional) Robust Analysis
if run_robust
    disp('---------------------------------------------');
    disp('Running robust per-animal analysis...');

    % Aggregate per animal
    animal_ID = lme_model.Variables.animal_ID; % Assumes these names
    group = lme_model.Variables.Group;
    response = lme_model.Variables.Response;

    [G, animal_groups] = findgroups(animal_ID);
    mean_response_per_animal = splitapply(@nanmean, response, G);
    genotype_per_animal = splitapply(@(x) x(1), group, G);

    % Design matrix
    X = double(genotype_per_animal == 'OxtKO');
    [b, stats] = robustfit(X, mean_response_per_animal);

    disp(['Robust Genotype Effect Estimate: ', num2str(b(2))]);
    disp(['Robust p-value for Genotype: ', num2str(stats.p(2))]);

    % Add robust results to summary text
    summary_text = [summary_text; {sprintf('Robust beta = %.3f', b(2))}];
    summary_text = [summary_text; {sprintf('Robust p = %.4g', stats.p(2))}];
end

%% Step 6: Plot summary annotation
annotation('textbox', [0.15 0.82 0.25 0.15], 'String', summary_text, ...
    'FitBoxToText', 'on', 'BackgroundColor', box_color, ...
    'EdgeColor', 'k', 'FontSize', 10, 'Interpreter', 'tex');

%% Step 7: Save Plot if requested
if ~isempty(plot_save_path)
    fprintf('Saving diagnostic plots to: %s\n', plot_save_path);
    exportgraphics(fig, fullfile(plot_save_path, ['LMEcheck_' fig_name '.pdf']), 'ContentType','vector');
end

%% Step 8: Console Output
disp('---------------------------------------------');
disp('Formal Residual Slope Test (Residuals ~ Fitted Values)');
disp(mdl_res.Coefficients);

if pval_slope < 0.05
    disp(['⚠️ WARNING: Significant residual trend detected (p = ', num2str(pval_slope), ').']);
    disp(['Slope estimate: ', num2str(slope_estimate)]);
    disp('Residuals depend systematically on fitted values.');
else
    disp(['✅ OK: No significant residual trend (p = ', num2str(pval_slope), ').']);
end
disp('---------------------------------------------');

fprintf('\nSummary of LME Model Assumption Checks:\n');
fprintf(' - Linearity: %s\n', ternary(pval_slope < 0.05, 'Potential violation (trend detected)', 'OK'));
fprintf(' - Homoscedasticity: Visually check Residuals vs Fitted (should be random scatter)\n');
fprintf(' - Normality: Visually check QQ Plot and Histogram\n');

end

%% Helper function
function out = ternary(condition, valTrue, valFalse)
    if condition
        out = valTrue;
    else
        out = valFalse;
    end
end
