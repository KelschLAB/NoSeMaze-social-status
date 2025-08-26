function generate_correlation_matrix(data_table, NoSeMaze_metrics, Basic_metrics, day_ranges, output_dir)
% Generate and save a correlation matrix between NoSeMaze_metrics and Basic_metrics.
%
% Parameters:
% data_table - Summary data table
% NoSeMaze_metrics - Cell array of NoSeMaze metric names
% Basic_metrics - Cell array of basic metric names
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
    num_NoSeMaze_metrics = numel(NoSeMaze_metrics);
    num_Basic_metrics = numel(Basic_metrics);
    corr_matrix = nan(num_NoSeMaze_metrics, num_Basic_metrics);
    pval_matrix = nan(num_NoSeMaze_metrics, num_Basic_metrics);
    
    % Compute correlations
    for i = 1:num_NoSeMaze_metrics
        for j = 1:num_Basic_metrics
            xData = current_data.(NoSeMaze_metrics{i});
            yData = current_data.(Basic_metrics{j});
            
            if isnumeric(yData) % Exclude non-numeric variables
                [rho, pval] = corr(xData, yData, 'Rows', 'complete', 'Type', 'Pearson');
            else
                rho = nan; % Correlation not applicable
                pval = nan;
            end
            
            corr_matrix(i, j) = rho;
            pval_matrix(i, j) = pval;
        end
    end
    
    % Save correlation matrix to a CSV file
    corr_table = array2table(corr_matrix, 'RowNames', NoSeMaze_metrics, 'VariableNames', Basic_metrics);
    writetable(corr_table, fullfile(output_dir, sprintf('CorrelationMatrix_%s.csv', range_name)), 'WriteRowNames', true);
    
    % Plot heatmap of correlation matrix
    fig = figure('Visible','off');
    imagesc(corr_matrix);
    % Axis
    ax=gca;
    axis square
    ax.XTick=[1:length(Basic_metrics)];
    ax.XTickLabel=Basic_metrics;
    ax.XTick=[1:length(NoSeMaze_metrics)];
    ax.YTickLabel=NoSeMaze_metrics;
    ax.TickLabelInterpreter = 'none';
    
    % title
    tt=title(sprintf('Correlation Matrix (%s)', range_name));
    tt.Interpreter='none';
    
    % colorbar
    ax.Colormap=redbluecmap;
    ax.CLim=[-1,1];
    colorbar;
    
    % Overlay correlation coefficients and significance markers
    hold on;
    for i = 1:num_NoSeMaze_metrics
        for j = 1:num_Basic_metrics
            if pval_matrix(i, j) < 0.05
                % Display both r and asterisk for significant values
                text(j, i, sprintf('%.2f*', corr_matrix(i, j)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold');
            else
                % Display r without asterisk for non-significant values
                text(j, i, sprintf('%.2f', corr_matrix(i, j)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Color', 'black', 'FontSize', 10, 'FontWeight', 'normal');
            end
        end
    end
    hold off
    
    % Save the heatmap
    [annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
    exportgraphics(fig, fullfile(output_dir, sprintf('CorrelationMatrix_%s.pdf', range_name)));
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
