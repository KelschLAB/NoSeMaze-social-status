function [data_rep1,data_rep2,animal_ID,genotype_animals,results_table]=plot_and_analyze_repetition_stability(data_table, metric, day_ranges, output_dir, LME_covariates, rep_first, rep_second)
% Analyze stability of a given metric across repetitions and day ranges
%
% Parameters:
% data_table - Summary data table
% metric - Metric to analyze (e.g., 'Rank_Competition')
% day_ranges - Cell array of day ranges to analyze (e.g., {'D1_14', 'D1_21'})
% output_dir - Directory to save results
% LME_covariates - Optional covariates for LME model

LME_covariates_names = strjoin(LME_covariates, '_');

%% Iterate Over Day Ranges
for day_idx = 1:numel(day_ranges)
    range_name = day_ranges{day_idx};

    % Extract relevant data for the current day range
    if ~isfield(data_table, range_name)
        warning('Day range %s not found in data_table. Skipping...', range_name);
        continue;
    end
    current_data = data_table.(range_name);

    % Get unique animals and their repetitions
    unique_animals = unique(current_data.Mouse_RFID);
    repetitions_index = cell(size(unique_animals));

    for an = 1:numel(unique_animals)
        animal_rows = find(strcmp(current_data.Mouse_RFID, unique_animals{an}));
        repetitions_index{an} = animal_rows;
    end

    %% Analyze Stability for Repetitions
    ii = rep_first; % First repetition
    jj = rep_second; % Second repetition
    % xData = [];
    % yData = [];
    all_data = nan(length(unique_animals),max(cellfun(@numel,repetitions_index)));
    %% Prepare genotype data (handles empty LME_covariates)
    if ~isempty(LME_covariates) && ~strcmp(LME_covariates{1}, '')
        genotype_field = LME_covariates{1};
        genotype_data = nan(length(unique_animals), max(cellfun(@numel,repetitions_index)));
        genotype_original = cell(length(unique_animals), max(cellfun(@numel,repetitions_index)));

        for an = 1:numel(unique_animals)
            for rep = 1:length(repetitions_index{an})
                genotype_data(an,rep) = strcmp(current_data.(genotype_field)(repetitions_index{an}(rep)),'OxtKO');
                genotype_original(an,rep) = current_data.(genotype_field)(repetitions_index{an}(rep));
            end
        end
    else
        % If no covariate is provided, fill with default zeros / empty cells
        genotype_data = zeros(length(unique_animals), max(cellfun(@numel,repetitions_index)));
        genotype_original = repmat({''}, length(unique_animals), max(cellfun(@numel,repetitions_index)));
    end

    counter = 1;

    for an = 1:numel(unique_animals)
        for rep = 1:length(repetitions_index{an})
            % Always assign metric
            all_data(an,rep) = current_data.(metric)(repetitions_index{an}(rep));

            % Only assign genotype/covariate if it exists
            if ~isempty(LME_covariates) && ~strcmp(LME_covariates{1}, '')
                genotype_data(an,rep) = strcmp(current_data.(LME_covariates{1})(repetitions_index{an}(rep)), 'OxtKO');
                genotype_original(an,rep) = current_data.(LME_covariates{1})(repetitions_index{an}(rep));
            else
                % Default values if no covariate is provided
                genotype_data(an,rep) = 0;
                genotype_original(an,rep) = {''};
            end
        end
        selection_idx(an) = repetitions_index{an}(ii); % keep this for LME indexing
    end


    % Apply function to move non-NaN values to the beginning of each row,
    % while keeping the order of appearance, and pushing NaNs to the end.
    sortedMatrix = sortNonNaNFirst(all_data);

    % Data
    % valid_idx = ~isnan(sortedMatrix(:,ii)) & ~isnan(sortedMatrix(:,jj));
    % xData = sortedMatrix(valid_idx,ii);
    % yData = sortedMatrix(valid_idx,jj);
    % animal_ID = unique_animals(valid_idx);
    % selection_idx = selection_idx(valid_idx);
    xData = sortedMatrix(:,ii);
    yData = sortedMatrix(:,jj);
    animal_ID = unique_animals(:);
    selection_idx = selection_idx(:);
    genotype_animals = genotype_original(:,1);

    % Initialize LME Table
    if ~isempty(LME_covariates) && ~strcmp(LME_covariates{1}, '')
        % Initialize with Animal ID, xData, and yData
        LME_DataTable = table(animal_ID, xData, yData, 'VariableNames', {'Mouse_RFID', [metric '_R1'], [metric '_R2']});

        % Add covariates dynamically
        for lx = 1:length(LME_covariates)
            covariate_name = LME_covariates{lx};
            LME_DataTable.(covariate_name) = current_data.(covariate_name)(selection_idx);
        end
    else
        % Only include Animal ID, xData, and yData if no covariates
        LME_DataTable = table(animal_ID, xData, yData, 'VariableNames', {'Mouse_RFID', [metric '_R1'], [metric '_R2']});
    end

    % Count the number of entries in each cell
    numEntries = cellfun(@(x) numel(x), repetitions_index);

    % Find included animals
    included_animals = unique_animals(find(numEntries >= 2));

    %% Plot Correlation
    fig = figure('visible','off');
    %         set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.9]);
    % Correlation coefficients
    [rho, pval] = corr(xData, yData, 'type', 'Pearson', 'Rows', 'pairwise');
    [rho_sp, pval_sp] = corr(xData, yData, 'type', 'Spearman', 'Rows', 'pairwise');
    [rho_partial, pval_partial] = partialcorr(xData, yData,  genotype_data(:,1), 'type', 'Pearson', 'Rows', 'pairwise');
    [rho_sp_partial, pval_sp_partial] = partialcorr(xData, yData, genotype_data(:,1), 'type', 'Spearman', 'Rows', 'pairwise');
    % Scatter plot
    sc=scatter(xData, yData);
    sc.SizeData=20;
    sc.MarkerEdgeColor='none';
    sc.MarkerFaceColor=[0,0,0];
    hold on;
    % Customize axis
    box('off');
    ax=gca;
    axis square;
    % correlation line
    ll=lsline;
    ll.Color=[0,0,0];
    ll.LineWidth=2;

    if contains(metric,'Rank') || contains(metric,'N_')
        ax.XLim(1)=0;
        ax.YLim(1)=0;
    elseif contains(metric,'cuberoot_Fraction')
        ax.YLim=[0,1];
        ax.XLim=[0,1];
    elseif contains(metric,'cuberoot_') || min([min(xData),min(yData)])>=0
        max_all_data = max([abs(min(xData)),abs(min(yData)),max(xData),max(yData)]);
        ax.XLim=([0 1.1*max_all_data]);
        ax.YLim=([0 1.1*max_all_data]);
    else
        max_all_data = max([abs(min(xData)),abs(min(yData)),max(xData),max(yData)]);
        ax.XLim=([-1.1*max_all_data 1.1*max_all_data]);
        ax.YLim=([-1.1*max_all_data 1.1*max_all_data]);
    end
    ax.LineWidth=1.5;
    ax.FontSize=16;


    % Highlight overlapping points
    A = [xData,yData];
    A = A(~any(isnan(A),2),:);
    un_combs = unique(A,'rows');
    for ab = 1:size(un_combs,1)
        dot_size(ab) = 30*(nnz(sum(A==un_combs(ab,:),2)==2)-1)+3;
        scatter(un_combs(ab,1),un_combs(ab,2),dot_size(ab),'k','filled');
    end

    % Customize Plot
    ax.XLabel.String = {metric,['round ' num2str(ii)]};
    ax.XLabel.Interpreter= 'none';
    ax.YLabel.String = {metric,['round ' num2str(jj)]};
    ax.YLabel.Interpreter= 'none';
    title({sprintf('Stability: %s (%s)', metric, range_name), ...
        sprintf('Pearson R = %.2f, p = %.3f', rho, pval), ...
        sprintf('Spearman R = %.2f, p = %.3f', rho_sp, pval_sp)},'Interpreter','none');

    %% Compute statistics
    [~, p_values, t_sums, ~ ] = permutest( xData(~isnan(xData) & ~isnan(yData))', yData(~isnan(xData) & ~isnan(yData))', true, ...
        0.05, 10000, true);
    if isempty(t_sums); t_sums = NaN; end
    [p_wilcoxon, ~, stats] = signrank(xData(~isnan(xData) & ~isnan(yData))', yData(~isnan(xData) & ~isnan(yData))');
    [p_value_med, diff_obs_med] = paired_permutation_test(xData(~isnan(xData) & ~isnan(yData))', yData(~isnan(xData) & ~isnan(yData))', 10000, 'median', 'two-tailed');
    disp([metric ': signedrank=' num2str(stats.signedrank) ', p=' num2str(p_wilcoxon) ', diffPERM=' num2str(diff_obs_med) ', pPERM=' num2str(p_value_med)]);

    % Create a table with results
    % Ensure all values are in a **single row (1xN format)**
    values = [rho_sp, pval_sp, rho, pval, t_sums, p_values, stats.signedrank, p_wilcoxon, diff_obs_med, p_value_med];

    % Define column names (metric names)
    column_names = {'r_spearman', 'pval_spearman', 'r_pearson', 'pval_pearson', ...
        't_perm_cluster', 'pval_perm_cluster', 'z_wilcoxon', 'pval_wilcoxon', ...
        'diff_obs_median', 'pval_perm_median'};

    % Create the table with **one row** and multiple **column variables**
    results_table.(range_name) = array2table(values, 'VariableNames', column_names, 'RowNames', string(metric));

    %% Display statistical results in the figure
    stats_text = {
        [metric ' Results:']
        ['signedrank = ' num2str(stats.signedrank)]
        ['p = ' num2str(p_wilcoxon)]
        ['diffPERM = ' num2str(diff_obs_med)]
        ['pPERM (Median) = ' num2str(p_value_med)]
        ['rho_sp_partial = ' num2str(rho_sp_partial)]
        ['p_sp_partial = ' num2str(pval_sp_partial)]
        };

    % Position the text in the upper-right of the figure
    xlims = xlim;
    ylims = ylim;
    text(xlims(2) * 0.8, ylims(2) * 0.9, stats_text, 'FontSize', 8, 'BackgroundColor', 'w', 'Interpreter', 'none');

    %% Save Plot
    [annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
    exportgraphics(fig, fullfile(output_dir, sprintf('%s_correlation_repetition_%d_to_%d_cov_%s_%s.pdf', metric, ii, jj, LME_covariates_names, range_name)), 'Resolution', 300);

    %% Save Source Data
    source_data = table(xData, yData, ...
        'VariableNames', {sprintf('%s_R%d_%s', metric, ii, range_name), ...
        sprintf('%s_R%d_%s', metric, jj, range_name)});
    % Save source data
    writetable(source_data, fullfile(output_dir, sprintf('SourceData_%s_correlation_repetition_%d_to_%d_cov_%s_%s.csv', metric, ii, jj, LME_covariates_names, range_name)));

    %% LME
    LME_DataTable.Mouse_RFID = categorical(LME_DataTable.Mouse_RFID);

    % Define the dependent variable and grouping variable
    dependent_variable = [metric '_R2']; % Target variable for the LME
    grouping_variable = 'Mouse_RFID';             % Random effects grouping variable

    % Get all variable names from the table
    all_variables = LME_DataTable.Properties.VariableNames;

    % Exclude the dependent variable and grouping variable from the predictors
    predictors = setdiff(all_variables, {dependent_variable, grouping_variable});

    % Create the fixed-effects part of the formula
    fixed_effects = strjoin(predictors, ' + '); % Combine all predictors with ' + '

    % Build the complete formula
    %             formula = sprintf('%s ~ %s + (%s|%s)', dependent_variable, fixed_effects, predictors{1}, grouping_variable);
    formula = sprintf('%s ~ %s + (1|%s)', dependent_variable, fixed_effects, grouping_variable);

    % Display the formula (for verification)
    disp(['Constructed LME formula: ', formula]);

    % Fit the linear mixed-effects model
    lme = fitlme(LME_DataTable, formula);

    % Save SourceData and LME fixed effects results to a CSV file
    fixedEffectsTable = lme.Coefficients;
    fixedEffectsTable = table(lme.CoefficientNames',fixedEffectsTable.Estimate, fixedEffectsTable.SE, fixedEffectsTable.tStat, fixedEffectsTable.pValue, fixedEffectsTable.Lower, fixedEffectsTable.Upper, ...
        'VariableNames', {'Names','Estimate', 'SE', 'tStat', 'pValue', 'CI_lower', 'CI_upper'});
    writetable(fixedEffectsTable, fullfile(output_dir, sprintf('LME_fixedEffects_%s_%d_to_%d_cov_%s_%s.csv', metric, ii, jj, LME_covariates_names, range_name)));
    % Save source data
    writetable(fixedEffectsTable, fullfile(output_dir, sprintf('LME_fixedEffects_%s_%d_to_%d_cov_%s_%s.csv', metric, ii, jj, LME_covariates_names, range_name)));



    % Save LME Results
    writetable(LME_DataTable, fullfile(output_dir, sprintf('SourceData_%s_LME_repetition_%d_to_%d_cov_%s_%s.csv', metric, ii, jj, LME_covariates_names, range_name)));

    close(fig);

    % define output
    data_rep1.(range_name) = xData;
    data_rep2.(range_name) = yData;
end
end

function [p_value, diff_obs] = paired_permutation_test(x, y, n_perm, method, tail)
    % Ensures x and y are column vectors
    x = x(:);
    y = y(:);

    % Check method input
    if strcmp(method, 'mean')
        diff_obs = mean(x - y); % Observed mean difference
        compute_stat = @mean;   % Function handle for permutation
    elseif strcmp(method, 'median')
        diff_obs = median(x - y); % Observed median difference
        compute_stat = @median;   % Function handle for permutation
    else
        error('Invalid method. Choose either ''mean'' or ''median''.');
    end

    % Default tail option
    if nargin < 5
        tail = 'two-tailed'; % Default to two-tailed test
    end

    % Permutation test
    n = length(x);
    perm_diffs = zeros(n_perm, 1);

    for i = 1:n_perm
        flip_signs = (rand(n, 1) > 0.5) * 2 - 1; % Randomly flip signs
        perm_diffs(i) = compute_stat((x - y) .* flip_signs);
    end

    % Compute p-value based on tail option
    switch tail
        case 'two-tailed'
            p_value = mean(abs(perm_diffs) >= abs(diff_obs)); % Two-tailed test
        case 'greater'
            p_value = mean(perm_diffs >= diff_obs); % One-tailed (x > y)
        case 'less'
            p_value = mean(perm_diffs <= diff_obs); % One-tailed (x < y)
        otherwise
            error('Invalid tail option. Choose ''two-tailed'', ''greater'', or ''less''.');
    end
end

function sortedMatrix = sortNonNaNFirst(matrix)
    % Function to move non-NaN values to the beginning of each row,
    % while keeping the order of appearance, and pushing NaNs to the end.
    
    [numRows, numCols] = size(matrix);  % Get matrix size
    sortedMatrix = NaN(numRows, numCols);  % Initialize with NaNs

    for i = 1:numRows
        rowValues = matrix(i, :);  % Extract row
        nonNaNValues = rowValues(~isnan(rowValues));  % Extract non-NaN values
        
        % Place non-NaN values at the beginning, keep NaNs at the end
        sortedMatrix(i, 1:length(nonNaNValues)) = nonNaNValues;
    end
end
