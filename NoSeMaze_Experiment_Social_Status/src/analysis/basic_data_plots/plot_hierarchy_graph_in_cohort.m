function f = plot_hierarchy_graph_in_cohort(match_matrix, data, IDs, cohort, metric_name)

[win, los]=find(match_matrix);
for i=1:length(win)
    weight(i)=match_matrix(win(i),los(i));
end
[~, sort_idx] = sort(data, 'descend');
[~, sort_idx_2] = sort(sort_idx, 'ascend');
%%%%% Graph plot
G=digraph(win,los,weight);  
LWidths = 3*G.Edges.Weight/max(G.Edges.Weight);
f = figure('Name', metric_name,'Visible','off');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.3, 0.4]);

p=plot(G,'LineWidth',LWidths,'NodeLabel',[],'Layout','force');
p.ArrowSize=8; p.EdgeColor='b';p.NodeColor=[.6 .6 .6];
p.NodeLabel=IDs;
p.NodeFontName = 'Arial';
p.NodeFontSize = 12;
p.MarkerSize = 10;
% color the nodes with rank
p.NodeCData = sort_idx_2;
% p.NodeCData = 10:10:100;
% colormap('jet');
% Flip the colormap
colormap(flipud(colormap));  % This flips the current colormap
% Add a colorbar
c = colorbar;
c.Label.String = 'rank';
set(c, 'YDir', 'reverse');

set_fonts()
f.Units = 'centimeters';
f.Position = [3 3 12 6];
title(['GraphPlot_',metric_name,'_', cohort],'Interpreter','none')

end