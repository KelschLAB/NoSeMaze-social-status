function plot_and_analyze_metric(...
    data_table, metric, day_ranges, output_dir, LME_covariates)
    % Updated to handle consistent structures for correlation and LME results

    % Initialize structures
    corrResults = struct('range1', {}, 'range2', {}, 'r_p', {}, 'p_p', {});
    lmeResults = struct('range1', {}, 'range2', {}, 'Formula', {}, 'FixedEffects', {}, 'RandomEffects', {}, 'ModelFit', {});

    % Number of comparisons
    comparison_count = (length(day_ranges) * (length(day_ranges) - 1)) / 2;
    fig = figure('Visible', 'on');
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.9]);
    
    subpl_counter = 1;

    for range_idx_1 = 1:length(day_ranges)
        for range_idx_2 = range_idx_1+1:length(day_ranges)
            range_name_1 = day_ranges{range_idx_1};
            range_name_2 = day_ranges{range_idx_2};
            
            % Extract data for current day ranges
            xData = data_table.(range_name_1).(metric);
            yData = data_table.(range_name_2).(metric);
            
            %% Create subplot
            % scatter plot
            subplot(1, comparison_count, subpl_counter);
            sc=scatter(xData, yData);            
            sc.SizeData=20;
            sc.MarkerEdgeColor='none';
            sc.MarkerFaceColor=[0,0,0];
            
            % Customize axis
            box('off');
            ax = gca;
            axis square;
            ax.XLabel.String = range_name_1;
            ax.YLabel.String = range_name_2;
            ax.XLabel.Interpreter = 'none';
            ax.YLabel.Interpreter = 'none';
            ax.FontSize = 10;
            ax.LineWidth = 1;
            if contains(metric,'Rank') || contains(metric,'N_') 
                ax.YLim(1) = 0;
                ax.XLim(1) = 0;
            elseif contains(metric,'cuberoot_Fr')
                ax.YLim = [0,1];
                ax.XLim = [0,1];
            end
            hold on;
            
            % correlation line
            ll=lsline;
            ll.Color=[0,0,0];
            ll.LineWidth=1.5;
            
            % adapt dot-size for rank
            if contains(metric,'Rank') || contains(metric,'N_') || contains(metric,'cuberoot_Fr') || contains(metric,'boxcox_Fr')
                hold on
                A = [xData,yData];
                A = A(~any(isnan(A), 2), :);
                un_combs = unique(A,'rows');
                for ab = 1:size(un_combs,1)
                    dot_size(ab) = 30*(nnz(sum(A==un_combs(ab,:),2)==2)-1)+10;
                    scatter(un_combs(ab,1),un_combs(ab,2),dot_size(ab),'k','filled');
                end
            end
            
            % title
            title(sprintf('%s: %s vs %s', metric, range_name_1, range_name_2), 'Interpreter', 'none');
         
            % Calculate Pearson and Spearman correlation and annotate plot 
            [r_p, p_p] = corr(xData, yData, 'type', 'Pearson', 'rows', 'complete');
            [r_sp, p_sp] = corr(xData, yData, 'type', 'Spearman', 'rows', 'complete');
            tx1 = text(ax.XLim(1) + .1 * (diff(ax.XLim)), ax.YLim(1) + .9 * (diff(ax.YLim)), ['r_p = ' num2str(r_p) ', p_p = ' num2str(p_p)]);
            if p_p < .05
                tx1.Color = [0.8, 0, 0];
            end
            tx2 = text(ax.XLim(1) + .1 * (diff(ax.XLim)), ax.YLim(1) + .75 * (diff(ax.YLim)), ['r_s_p = ' num2str(r_sp) ', p_s_p = ' num2str(p_sp)]);
            if p_sp < .05
                tx2.Color = [0.8, 0, 0];
            end      
            
            % Save correlation SourceData and CorrelationResults to a CSV file
            sourceDataTable = table(data_table.(range_name_1).(metric), data_table.(range_name_2).(metric), data_table.(range_name_1).Mouse_RFID, data_table.(range_name_1).cohort, 'VariableNames', {range_name_1,range_name_2,'Mouse_RFID','cohort'});
            writetable(sourceDataTable, fullfile(output_dir, ['SourceData_correlation_' metric '_'  range_name_1 '_to_' range_name_2 '.csv']));
            
            correlationResults = table(r_p, p_p, r_sp, p_sp, 'VariableNames', {'r_Pearson', 'p_Pearson', 'r_Spearman', 'p_Spearman'});
            writetable(correlationResults, fullfile(output_dir, ['CorrelationResults_' metric '_'  range_name_1 '_to_' range_name_2 '_.csv']));
            
            %% LME
            % Ensure cohort is a categorical variable
            if ~isempty(LME_covariates) && ~strcmp(LME_covariates{1}, '')
                LME_DataTable_1 = data_table.(range_name_1)(:, [metric, 'Mouse_RFID', 'cohort', LME_covariates]);
                LME_DataTable_1.Properties.VariableNames = [[metric '_' range_name_1], 'Mouse_RFID', 'cohort', LME_covariates];
                LME_covariates_names = strjoin(LME_covariates, '_');
            else
                LME_DataTable_1 = data_table.(range_name_1)(:, {metric, 'Mouse_RFID', 'cohort'});
                LME_DataTable_1.Properties.VariableNames = {[metric '_' range_name_1], 'Mouse_RFID', 'cohort'};
                LME_covariates_names = '';
            end
            LME_DataTable_2 = data_table.(range_name_2)(:, metric);
            LME_DataTable_2.Properties.VariableNames = {[metric '_' range_name_2]};
            LME_DataTable = [LME_DataTable_2,LME_DataTable_1];
            
            LME_DataTable.Mouse_RFID = categorical(LME_DataTable.Mouse_RFID);
            
            % Define the dependent variable and grouping variable
            dependent_variable = [metric '_' range_name_2]; % Target variable for the LME
            grouping_variable = 'Mouse_RFID';             % Random effects grouping variable
            
            % Get all variable names from the table
            all_variables = LME_DataTable.Properties.VariableNames;
            
            % Exclude the dependent variable and grouping variable from the predictors
            predictors = setdiff(all_variables, {dependent_variable, grouping_variable, 'cohort'});
            
            % Create the fixed-effects part of the formula
            fixed_effects = strjoin(predictors, ' + '); % Combine all predictors with ' + '
            
            % Build the complete formula
%             formula = sprintf('%s ~ %s + (%s|%s)', dependent_variable, fixed_effects, predictors{1}, grouping_variable);
            formula = sprintf('%s ~ %s + (1|%s)', dependent_variable, fixed_effects, grouping_variable);
            
            % Display the formula (for verification)
            disp(['Constructed LME formula: ', formula]);
            
            % Fit the linear mixed-effects model
            lme = fitlme(LME_DataTable, formula);
            
%             % Display results
%             disp(lme);
            
            % Save SourceData and LME fixed effects results to a CSV file
            fixedEffectsTable = lme.Coefficients;
            fixedEffectsTable = table(lme.CoefficientNames',fixedEffectsTable.Estimate, fixedEffectsTable.SE, fixedEffectsTable.tStat, fixedEffectsTable.pValue, fixedEffectsTable.Lower, fixedEffectsTable.Upper, ...
                'VariableNames', {'Names','Estimate', 'SE', 'tStat', 'pValue', 'CI_lower', 'CI_upper'});
            writetable(fixedEffectsTable, fullfile(output_dir, ['LMEfixedEffects_' metric '_'  range_name_1 '_to_' range_name_2 '_cov_' LME_covariates_names '.csv']));
            % Convert Mouse_RFID in LME_DataTable to cell array of char if needed
            if ~iscell(LME_DataTable.Mouse_RFID)
                LME_DataTable.Mouse_RFID = cellstr(string(LME_DataTable.Mouse_RFID));
            end
            writetable(join(data_table.(range_name_1)(:,{'Mouse_RFID','cohort'}),LME_DataTable,'Keys',{'Mouse_RFID','cohort'}), fullfile(output_dir, ['SourceData_LME_' metric '_'  range_name_1 '_to_' range_name_2 '_cov_' LME_covariates_names '.csv']));

            % update the counter
            subpl_counter = subpl_counter + 1;                
        end
    end
    
    % Save figure
    [annot, srcInfo] = docDataSrc(fig, fullfile(output_dir), mfilename('fullpath'), logical(1));
    exportgraphics(fig, fullfile(output_dir, sprintf('%s_scatterplots.pdf', metric)), 'ContentType','vector','Resolution', 300);
    close(fig);
end
