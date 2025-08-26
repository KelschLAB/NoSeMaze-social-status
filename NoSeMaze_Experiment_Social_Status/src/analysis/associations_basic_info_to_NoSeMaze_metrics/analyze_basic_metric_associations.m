function analyze_basic_metric_associations(data_table, NoSeMaze_metric, Basic_metrics, day_ranges, output_dir, LME_covariates)
% Analyze associations between Basic_metrics and a specified NoSeMaze_metric.
%
% Parameters:
% data_table - Summary data table
% NoSeMaze_metric - NoSeMaze_metric to analyze (e.g., 'DSz_Competition')
% Basic_metrics - Cell array of basic_metric names (e.g., {'weight_pre', 'genotype'})
% day_ranges - Day ranges to analyze (e.g., {'D1_14', 'D1_21'})
% output_dir - Directory to save results

for day_idx = 1:numel(day_ranges)
    range_name = day_ranges{day_idx};
    
    % Check if day range exists in the data table
    if ~isfield(data_table, range_name)
        warning('Day range %s not found in data table. Skipping...', range_name);
        continue;
    end
    
    % Extract data for the current day range
    current_data = data_table.(range_name);
    
    % Prepare data for analysis
    xData = current_data.(NoSeMaze_metric);
    % Initialize LME Table
    if ~isempty(LME_covariates) && ~strcmp(LME_covariates{1}, '')
        % Initialize with Animal ID, xData, and yData
        LME_DataTable = table(current_data.Mouse_RFID, xData, 'VariableNames', {'Mouse_RFID', NoSeMaze_metric});

        % Add covariates dynamically
        for lx = 1:length(LME_covariates)
            covariate_name = LME_covariates{lx};
            LME_DataTable.(covariate_name) = current_data.(covariate_name);
        end
    else
        % Only include Animal ID, xData, and yData if no covariates
        LME_DataTable = table(current_data.Mouse_RFID, xData, 'VariableNames', {'Mouse_RFID', NoSeMaze_metric});
    end
    
    % Add Basic_metrics to the table
    for basic_metric_idx = 1:numel(Basic_metrics)
        basic_metric_name = Basic_metrics{basic_metric_idx};
        LME_DataTable.(basic_metric_name) = current_data.(basic_metric_name);
    end
    
    % Display number of entries for debugging
    fprintf('Day range: %s | NoSeMaze_metric: %s | Entries: %d\n', range_name, NoSeMaze_metric, height(LME_DataTable));
    
    %% Generate Scatter Plots and notBoxPlots
    for basic_metric_idx = 1:numel(Basic_metrics)
        
        basic_metric_name = Basic_metrics{basic_metric_idx};
        
        xData = LME_DataTable.(basic_metric_name);
        yData = LME_DataTable.(NoSeMaze_metric);
        
        %% Plot Correlation
        fig = figure('visible','off');
        
        % Correlation coefficients
        if length(unique(LME_DataTable.(basic_metric_name)))~=2
            [rho, pval] = corr(xData, yData, 'type', 'Pearson');
            [rho_sp, pval_sp] = corr(xData, yData, 'type', 'Spearman');
        end
        % Plot scatter plot or notBoxPlot for binary vs continuous basic_metrics
        if length(unique(LME_DataTable.(basic_metric_name)))==2
            % notBoxPlot            
            nb = notBoxPlot(LME_DataTable.(NoSeMaze_metric), categorical(LME_DataTable.(basic_metric_name)));
            % color definitions
            for ib = 1:length(nb)
                nb(ib).data.MarkerSize=8;
                nb(ib).data.MarkerEdgeColor='none';
                nb(ib).semPtch.EdgeColor='none';
                nb(ib).sdPtch.EdgeColor='none';
            end
            nb(1).data.MarkerFaceColor= [204/255 51/255 204/255];
            nb(1).mu.Color= [204/255 51/255 204/255];
            nb(1).semPtch.FaceColor= [221/255 118/255 221/255];
            nb(1).sdPtch.FaceColor= [247/255 221/255 247/255];
            % color definitions
            nb(2).data.MarkerFaceColor= [0 160/255 227/255];
            nb(2).mu.Color= [0 160/255 227/255];
            nb(2).semPtch.FaceColor= [90/255 194/255 237/255];
            nb(2).sdPtch.FaceColor= [211/255 239/255 250/255];
            % Customize axis
            box('off');   
            ax=gca;
            axis square;
        else
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
        end
        
        if length(unique(LME_DataTable.(basic_metric_name)))~=2
            if contains(NoSeMaze_metric,'Rank') || contains(NoSeMaze_metric,'N_')
                ax.YLim(1)=0;
                if min(xData)<0
                    ax.XLim=([1.1*min(xData) 1.1*max(xData)]);
                else
                    ax.XLim=([0.9*min(xData) 1.1*max(xData)]);
                end
            elseif contains(NoSeMaze_metric,'Competition') || contains(NoSeMaze_metric,'Fraction_')
                ax.YLim = [0 0.8];
                if min(xData)<0
                    ax.XLim=([1.1*min(xData) 1.1*max(xData)]);
                else
                    ax.XLim=([0.9*min(xData) 1.1*max(xData)]);
                end
            else
                ax.YLim=([1.1*min(yData) 1.1*max(yData)]);
                if min(xData)<0
                    ax.XLim=([1.1*min(xData) 1.1*max(xData)]);
                else
                    ax.XLim=([0.9*min(xData) 1.1*max(xData)]);
                end
            end            
        end
        
        if length(unique(LME_DataTable.(basic_metric_name)))~=2
            % correlation line
            ll=lsline;
            ll.Color=[0,0,0];
            ll.LineWidth=1.5;
        end
        
        % Highlight overlapping points
        if length(unique(LME_DataTable.(basic_metric_name)))~=2            
            A = [xData,yData];
            un_combs = unique(A,'rows');
            for ab = 1:size(un_combs,1)
                dot_size(ab) = 30*(nnz(sum(A==un_combs(ab,:),2)==2)-1)+3;
                scatter(un_combs(ab,1),un_combs(ab,2),dot_size(ab),'k','filled');
            end
        end
        
        % Customize plot
        ax.XLabel.String = basic_metric_name;
        ax.XLabel.Interpreter = 'none';
        ax.YLabel.String = NoSeMaze_metric;
        ax.YLabel.Interpreter = 'none';
        if length(unique(LME_DataTable.(basic_metric_name)))~=2
            title({sprintf('%s vs %s (%s)', NoSeMaze_metric, basic_metric_name, range_name), ...
                sprintf('Pearson R = %.2f, p = %.3f', rho, pval), ...
                sprintf('Spearman R = %.2f, p = %.3f', rho_sp, pval_sp)},'Interpreter','none');
        else
            title({sprintf('%s vs %s (%s)', NoSeMaze_metric, basic_metric_name, range_name)},'Interpreter','none')
        end
        
        % Save Plot
        [annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
        if length(unique(LME_DataTable.(basic_metric_name)))==2
                exportgraphics(fig, fullfile(output_dir, sprintf('NotBoxPlot_%s_vs_%s_%s.pdf', NoSeMaze_metric, basic_metric_name, range_name)));
        else
                exportgraphics(fig, fullfile(output_dir, sprintf('Scatterplot_%s_vs_%s_%s.pdf', NoSeMaze_metric, basic_metric_name, range_name)));
        end
        close(fig);
        
        %% Linear Mixed-Effects Model
        LME_DataTable.Mouse_RFID = categorical(LME_DataTable.Mouse_RFID);

        % Define the dependent variable and grouping variable
        dependent_variable = NoSeMaze_metric; % Target variable for the LME
        grouping_variable = 'Mouse_RFID';             % Random effects grouping variable
    
        % Get all variable names from the table
        all_variables = LME_DataTable.Properties.VariableNames;
            
        % Exclude the dependent variable and grouping variable from the predictors
        if ~isempty(LME_covariates) && ~strcmp(LME_covariates{1}, '')
            predictors = [basic_metric_name,LME_covariates];
        else
            predictors = basic_metric_name;
        end
        LME_covariates_names = strjoin(LME_covariates, '_');

        % Create the fixed-effects part of the formula
        if ~isempty(LME_covariates) && ~strcmp(LME_covariates{1}, '')
            fixed_effects = strjoin(predictors, ' + '); % Combine all predictors with ' + '
        else
            fixed_effects = predictors;
        end
        % Build the complete formula
        %             formula = sprintf('%s ~ %s + (%s|%s)', dependent_variable, fixed_effects, predictors{1}, grouping_variable);
        formula = sprintf('%s ~ %s + (1|%s)', dependent_variable, fixed_effects, grouping_variable);
    
        % Display the formula (for verification)
        disp(['Constructed LME formula: ', formula]);
        
        %         % Define the fixed-effects formula
        % fixed_effects = X_metric_name; % Combine predictors
        % formula = sprintf('%s ~ %s + (1|Mouse_RFID)', Y_metric, fixed_effects);
        
        % Fit LME model
        try
            lme = fitlme(LME_DataTable, formula);
        catch ME
            warning('LME failed for %s in %s: %s', dependent_variable, range_name, ME.message);
            continue;
        end

        % Save LME results
        fixedEffectsTable = lme.Coefficients;
        fixedEffectsTable = table(lme.CoefficientNames',fixedEffectsTable.Estimate, fixedEffectsTable.SE, fixedEffectsTable.tStat, fixedEffectsTable.pValue, fixedEffectsTable.Lower, fixedEffectsTable.Upper, ...
        'VariableNames', {'Names','Estimate', 'SE', 'tStat', 'pValue', 'CI_lower', 'CI_upper'});
        % Permutation test:
        % % Get per-animal summaries (e.g., mean or median latency)
        % LME_DataTable = LME_DataTable(~isnan(LME_DataTable.(dependent_variable)),:);
        % [G, animal] = findgroups(LME_DataTable.Mouse_RFID);
        % animal_summary = splitapply(@nanmean, LME_DataTable.(dependent_variable), G);
        % group_labels = splitapply(@(x) x(1), LME_DataTable.genotype, G); % one label per animal
        % % Then run a simple permutation test on animal_summary grouped by group_labels
        % [observed_diff1, ~, p_value1] = permutation_test_unpaired(animal_summary(strcmp(group_labels,'WT')), animal_summary(strcmp(group_labels,'OxtKO')), 10000, 'mean')
        % [observed_diff2, ~, p_value2] = permutation_test_unpaired(animal_summary(strcmp(group_labels,'WT')), animal_summary(strcmp(group_labels,'OxtKO')), 10000, 'median')

        if ~isempty(LME_covariates) && ~strcmp(LME_covariates{1}, '')
            % writetable([fixedEffectsTable,table([observed_diff1;p_value1],[observed_diff2;p_value2],'VariableNames',{'perm_mean[diff,p]','perm_median[diff,p]'})], fullfile(output_dir, sprintf('LME_results_%s_vs_%s_cov_%s_%s.csv', NoSeMaze_metric, basic_metric_name, LME_covariates_names, range_name)));
            writetable(fixedEffectsTable,fullfile(output_dir, sprintf('LME_results_%s_vs_%s_cov_%s_%s.csv', NoSeMaze_metric, basic_metric_name, LME_covariates_names, range_name)));
            % Save source data
            writetable(LME_DataTable(:,[grouping_variable, dependent_variable, predictors]), fullfile(output_dir, sprintf('SourceData_%s_vs_%s_cov_%s_%s.csv', NoSeMaze_metric, basic_metric_name, LME_covariates_names, range_name)));
        else
            writetable(fixedEffectsTable,fullfile(output_dir, sprintf('LME_results_%s_vs_%s_%s.csv', NoSeMaze_metric, basic_metric_name, range_name)));
            % Save source data
            writetable(LME_DataTable(:,[{grouping_variable}, {dependent_variable}, predictors(:)']), fullfile(output_dir, sprintf('SourceData_%s_vs_%s_%s.csv', NoSeMaze_metric, basic_metric_name, range_name)));
        end
        % disp(dependent_variable);
        % [fixedEffectsTable,table([observed_diff1;p_value1],[observed_diff2;p_value2],'VariableNames',{'perm_mean[diff,p]','perm_median[diff,p]'})]
    end
end

