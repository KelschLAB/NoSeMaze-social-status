function generate_and_plot_association_matrix_basic_metrics(data_table, X_metrics, Y_metrics, day_ranges, output_dir, method)
% Generate and save an association matrix between X_metrics and Y_metrics.
%
% Parameters:
% data_table - Summary data table
% X_metrics - Cell array of NoSeMaze metric names
% Y_metrics - Cell array of basic metric names
% day_ranges - Cell array of day ranges (e.g., {'D1_14', 'D1_21'})
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
    
    % Initialize correlation matrices
    num_X_metrics = numel(X_metrics);
    num_Y_metrics = numel(Y_metrics);
    corr_matrix = nan(num_X_metrics, num_Y_metrics);
    pval_matrix = nan(num_X_metrics, num_Y_metrics);
    
    % Compute correlations
    for i = 1:num_X_metrics
        for j = 1:num_Y_metrics
            xData = current_data.(X_metrics{i});
            yData = current_data.(Y_metrics{j});
            
            if isnumeric(xData) % Exclude non-numeric variables
                [rho, pval] = corr(xData, yData, 'Rows', 'complete', 'Type', method);
            elseif length(unique(xData))==2
                xData_values = unique(xData);
                [rho, pval] = corr(strcmp(xData,xData_values{1}), yData, 'Rows', 'complete', 'Type', method);
            else
                rho = nan; % Correlation not applicable
                pval = nan;
            end
            
            corr_matrix(i, j) = rho;
            pval_matrix(i, j) = pval;
        end
    end
    
    % Save correlation matrix to a CSV file
    corr_table = array2table(corr_matrix, 'RowNames', X_metrics, 'VariableNames', Y_metrics);    
    p_table = array2table(pval_matrix, 'RowNames', X_metrics, 'VariableNames', Y_metrics);
    writetable(corr_table, fullfile(output_dir, sprintf('SourceData_AssociationMatrix_CorrCoeff_%s.csv', range_name)), 'WriteRowNames', true);
    writetable(p_table, fullfile(output_dir, sprintf('SourceData_AssociationMatrix_PVal_%s.csv', range_name)), 'WriteRowNames', true);
    
    % Plot heatmap of Association matrix
    fig = figure('Visible','off');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.9, 0.9]);
    plot_matrix = corr_matrix;
    % plot_matrix(triu(true(size(plot_matrix)))) = 0;
    imagesc(plot_matrix);
    % Axis
    ax=gca;
    box off
    axis square
    ax.XTick=[1:length(Y_metrics)];
    ax.XTickLabel=Y_metrics;
    ax.YTick=[1:length(X_metrics)];
    ax.YTickLabel=X_metrics;
    ax.TickLabelInterpreter = 'none';
    
    % title
    tt=title(sprintf('Association Matrix (%s, %s)', range_name, method));
    tt.Interpreter='none';
    
    % colorbar
    crameri('vik'); %ax.Colormap=redbluecmap;
    ax.CLim=[-1,1];
    colorbar;
    % Replace upper triangular half (excluding the diagonal) with NaN
    % pval_matrix(triu(true(size(pval_matrix)))) = NaN;
    [pID,pN,qvalues] = FDR(pval_matrix(~isnan(pval_matrix)),0.05);
    pBonferroni = 0.05/sum(sum(~isnan(pval_matrix)));
    % % % for r1=1:size(plot_matrix)
    % % %     for r2=1:size(plot_matrix)
    % % %         if pval_matrix(r1,r2)<=pBonferroni
    % % %             tx=text(r1,r2,'*');
    % % %             tx.FontSize=20;
    % % %         end
    % % %     end
    % % % end>

    % Overlay Association coefficients and significance markers
    hold on;
    for i = 1:num_X_metrics
        for j = 1:num_Y_metrics
            if pval_matrix(i, j) < pBonferroni
                % Display both r and asterisk for significant values
                text(j, i, sprintf('%.2f*', plot_matrix(i, j)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Color', 'white', 'FontSize', 8, 'FontWeight', 'bold');
            else
                % Display r without asterisk for non-significant values
                text(j, i, sprintf('%.2f', plot_matrix(i, j)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Color', 'black', 'FontSize', 4, 'FontWeight', 'normal');
            end
        end
    end
    hold off
    
    % Save the heatmap
    [annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
    exportgraphics(fig, fullfile(output_dir, sprintf('AssociationMatrix_%s.pdf', range_name)));
    close(fig);
end
end


function cmap = redbluecmap(n)
% Generates a red-to-white-to-blue colormap.
% Input: n - Number of colors (optional; defaults to 256).
if nargin < 1, n = 256; end
cmap = [linspace(0, 1, n/2)', linspace(0, 1, n/2)', ones(n/2, 1); ...
    ones(n/2, 1), linspace(1, 0, n/2)', linspace(1, 0, n/2)'];
end
