function [f,sorted_IDs,sorted_data,sorted_ranks] = plot_hierarchy_scores_in_cohort(data, IDs, cohort, metric_name)

[sorted_data, sort_idx] = sort(data, 'descend');
sorted_ranks = sort(sort_idx,'ascend'); % Sort ranks based on sorted data

% Create the figure
f = figure('Name', metric_name,'Visible','off');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.3, 0.4]);


% Plot David's score on the left y-axis
yyaxis left
dp = plot(sorted_data, 'k--d', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.2 0.2 0.2]);
xlim([0.5, length(sorted_data) + 0.5]);
if contains(metric_name, 'DSz')
    ylim([min(sorted_data) - 1, max(sorted_data) + 1])
else
    ylim([floor(min(sorted_data)), ceil(max(sorted_data)*10)/10])
end
ylabel(metric_name,'Interpreter','none', 'FontSize', 16)
box off;

% Plot hierarchy ranks on the right y-axis
yyaxis right
% plot(sorted_ranks(end:-1:1), 'r-o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r', ...
%     'MarkerFaceColor', [0.9 0.2 0.2]);
% Add a line through the dots
hold on;
plot(1:length(sorted_ranks), sorted_ranks(end:-1:1), '-', 'LineWidth', 1.5, 'Color', [0,0,1,0.5]);  % Add line with black color
% Plot the sorted ranks with a color map
cmap = parula(length(sorted_ranks));  % Use 'jet' or any other colormap of choice
scatter(1:length(sorted_ranks), sorted_ranks(end:-1:1), 60, sorted_ranks(end:-1:1), 'filled');
% colormap(flipud(cmap));  % Flip the colormap to reverse the color direction
colorbar;  % Add colorbar
caxis([min(sorted_ranks), max(sorted_ranks)]);  % Adjust the color scaling
ylim([min(sorted_ranks) - 1, max(sorted_ranks) + 1])
yticks(1:length(sorted_data));
yticklabels([length(sorted_data):-1:1]);
ylabel({'Rank',['(based on ', metric_name, ')']}, 'FontSize', 16, 'Interpreter', 'none')

% Add x-axis ticks and labels
xticks(1:length(sorted_data));
sorted_IDs = IDs(sort_idx);
xticklabels(sorted_IDs);
xtickangle(45);
set(gca, 'FontSize', 14) % Font size for tick labels

% Set general figure properties
%     set_fonts()
f.Units = 'centimeters';
f.Position = [3, 3, 9, 7]; % Adjust size to fit both plots nicely
title([metric_name ' and Ranks in ', cohort],'Interpreter','none')

% Add a legend for clarity
legend({metric_name, 'Ranks'}, 'Location', 'NorthEast','Interpreter','none');
end
