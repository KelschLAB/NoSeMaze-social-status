%% cluster_lickport_data.m
% Author: Jonathan Reinwald
% Date: March 2025

% Pre-clearing
clear all; clc; close all;

% Set pathes
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory
output_dir = fullfile(main_dir,'results','figures','cross_cohort','clustering_lickport_data');
if ~isfolder(output_dir)
    mkdir(output_dir);
end
addpath(genpath(fullfile(main_dir,'analysis','scripts','clustering_lickport_data')));
addpath(genpath(fullfile(main_dir,'src','helpers')));

% Load summary data
load(fullfile(main_dir,'data','processed','cross_cohort_files','summary_data.mat'));

% Load cohort information
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));
% Select cohorts that are included in both tube and lickport analyses
cohortsTbl = cohortsTbl((cohortsTbl.use_lickport == 1), :);

%% 📊 Compute Mean Switch Latencies and add them to the summary_data
% Compute the mean switch latencies for CS+ and CS-data = summary_data.lickport.D1_End;
cs_minus_switch = nanmean([summary_data.lickport.D1_End.cs_minus_switch_latency_at_cs_rev1, ...
    summary_data.lickport.D1_End.cs_minus_switch_latency_at_cs_rev2, ...
    summary_data.lickport.D1_End.cs_minus_switch_latency_at_cs_rev3, ...
    summary_data.lickport.D1_End.cs_minus_switch_latency_at_cs_rev4], 2);

cs_plus_switch = nanmean([summary_data.lickport.D1_End.cs_plus_switch_latency_at_cs_rev1, ...
    summary_data.lickport.D1_End.cs_plus_switch_latency_at_cs_rev2, ...
    summary_data.lickport.D1_End.cs_plus_switch_latency_at_cs_rev3, ...
    summary_data.lickport.D1_End.cs_plus_switch_latency_at_cs_rev4], 2);
% Add to summary data
summary_data.lickport.D1_End.cs_plus_switch_mean = cs_plus_switch;
summary_data.lickport.D1_End.cs_minus_switch_mean = cs_minus_switch;
summary_data.lickport.D1_End.boxcox_cs_plus_switch_mean = boxcox(cs_plus_switch);
summary_data.lickport.D1_End.boxcox_cs_minus_switch_mean = boxcox(cs_minus_switch);
summary_data.lickport.D1_End.cuberoot_cs_plus_switch_mean = cs_plus_switch.^(1/3);
summary_data.lickport.D1_End.cuberoot_cs_minus_switch_mean = cs_minus_switch.^(1/3);
% save summary data
save(fullfile(main_dir,'data','processed','cross_cohort_files','summary_data.mat'),'summary_data');

%% Filter input data (cohort selection)
fields = fieldnames(summary_data.lickport);
for field_idx = 1:numel(fields)
    summary_data_reduced.lickport.(fields{field_idx}) = summary_data.lickport.(fields{field_idx})(ismember(summary_data.lickport.(fields{field_idx}).cohort,cohortsTbl.cohort),:);
end
data = summary_data_reduced.lickport.D1_End;

%% 📌 Extract Features for Clustering
cs_plus_switch = data.cs_plus_switch_latency_at_cs_median;%data.cs_plus_switch_latency_at_cs_median;%cs_plus_time_to_criterion_median;
cs_minus_switch = data.cs_minus_switch_latency_at_cs_median;%data.cs_minus_switch_latency_at_cs_median;%cs_minus_time_to_criterion_median;
pre_cs_licking = data.baseline_rate_mean_omitfirst;  % Pre-CS licking (impulsivity)
correct_rejections = data.correct_rejection_rate_allTrials; % Timeout measure

cs_plus_modulation = data.cs_plus_modulation_peak;
cs_minus_modulation = data.cs_minus_modulation_min;
correct_hits = data.correct_hit_rate_allTrials;

id = data.Mouse_RFID;
cohort = data.cohort;

% Combine into a feature matrix
features = [cs_plus_switch, cs_minus_switch, pre_cs_licking, correct_rejections, correct_hits, cs_plus_modulation];
feature_names = {'cs_plus_switch', 'cs_minus_switch', 'pre_cs_licking', 'correct_rejections', 'correct_hits',  'cs_plus_modulation'};

% Remove rows with NaN values
valid_idx = ~any(isnan(features), 2);
features = features(valid_idx, :);
id = id(valid_idx);
cohort = cohort(valid_idx);

%% 🔄 Standardization (Min-Max Normalization)
% Min-Max is the best option for k-mean clustering
% Note that k-mean clustering does not require normal distribution.
% Measuring its performance with the silhouette values demonstrated better
% performance when not tranforming the data into normal distribution.
features_norm = (features - min(features)) ./ (max(features) - min(features));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 📌 Determine Optimal Number of Clusters (Elbow Method)
max_clusters = 6; % Set max number of clusters to test
sse = zeros(1, max_clusters);
for k = 1:max_clusters
    [~, ~, sumd] = kmeans(features_norm, k, 'Replicates', 100, 'Display', 'off');
    sse(k) = sum(sumd);
end

% Plot elbow curve
fig = figure('Visible','on');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.7]);
plot(1:max_clusters, sse, '-o');
xlabel('Number of Clusters (k)');
ylabel('Sum of Squared Errors (SSE)');
title('Elbow Method for Optimal k');
grid on;
% Axis
ax=gca;
box off
% Source Data
SourceDataTable = array2table([[1:max_clusters]', sse'],'VariableNames',{'cluster number','SSE'});
writetable(SourceDataTable,fullfile(output_dir, 'SourceData_EllbowPlot.csv'));
% Save the heatmap
[annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(output_dir, sprintf('Ellbow_Plot_Kmeans.pdf')));
close(fig);

%% 📌 K-Means Clustering (Using Optimal k)
% k = 2 is the optimal value in our data
num_clusters = 3; % Change based on elbow method result
[idx_kmeans, centroids] = kmeans(features_norm, num_clusters, 'Replicates', 1000);

% Evaluate clustering with silhouette score
silhouette_kmeans = mean(silhouette(features_norm, idx_kmeans));

% %% 📌 Gaussian Mixture Model (GMM) Clustering
% gmm = fitgmdist(features_norm, num_clusters, 'Replicates', 1000);  % Fit GMM
% idx_gmm = cluster(gmm, features_norm);  % Assign clusters
% 
% % Extract posterior probabilities
% gm_probabilities = posterior(gmm, features_norm);  % Probability of belonging to each cluster
% 
% % Evaluate clustering with silhouette score
% silhouette_gmm = mean(silhouette(features_norm, idx_gmm));

%% 📌 K-Medioids Clustering (Using Optimal k)
% k = 2 is the optimal value in our data
[idx_kmedoids, medoids] = kmedoids(features_norm, num_clusters, 'Replicates', 1000, 'Distance', 'euclidean');

% Evaluate clustering with silhouette score
silhouette_kmedoids = mean(silhouette(features_norm, idx_kmedoids));

% %% 📌 DBSCAN Clustering (Using Optimal k)
% % Set parameters for DBSCAN
% epsilon = 3; % Adjust based on density (use k-distance plot to find best value)
% minPts = 5;    % Minimum points required to form a cluster
%
% % Perform DBSCAN clustering
% idx_dbscan = dbscan(features_norm, epsilon, minPts);

%% 📌 PCA for Visualization
[coeff, score, latent, tsquared, explained, mu] = pca(features_norm);
reduced_data = score(:, 1:3); % Take first two principal components

% Visualize Clustering
% Define Custom RGB Colors
cluster_colors = [113 189 134; 178 170 209; 218 147 74] / 255;  % Normalize to [0,1]
cluster_colors = [113 189 134; 178 170 209; 218 147 74; 255 0 0; 0 0 255] / 255;
% K-Means Cluster Plot with Custom Colors
fig = figure('Visible','on');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.6], 'Renderer', 'painters');
hold on;
for i = 1:num_clusters
    scatter3(reduced_data(idx_kmeans == i, 1), reduced_data(idx_kmeans == i, 2), reduced_data(idx_kmeans == i, 3),...
        90, cluster_colors(i, :), 'filled');
    % scatter3(reduced_data(idx_kmeans == i, 1), reduced_data(idx_kmeans == i, 2), reduced_data(idx_kmeans == i, 3),...
    % 50, cluster_colors(i, :), 'filled');
end
view(-45,45);
hold off;
title('K-Means Clustering (PCA-reduced Data)');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
grid on;
% Axis
ax=gca;
ax.FontSize=20;
% ax.LineWidth=1;
box off
% axis square;
% Source Data
SourceDataTable = array2table([idx_kmeans,reduced_data],'VariableNames',{'cluster idx','PC1','PC2','PC3'});
writetable(SourceDataTable,fullfile(output_dir, 'SourceData_KmeanClustering_PCA.csv'));
% Save the figure
[annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(output_dir, ['KmeanClustering_PC_visualization_clusters_' num2str(num_clusters) '_features_' num2str(length(feature_names)) '.pdf']),'ContentType','vector');;
close(fig);

% K-Medoids Cluster Plot with Custom Colors
fig = figure('Visible','on');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.4]);
hold on;
for i = 1:num_clusters
    scatter3(reduced_data(idx_kmedoids == i, 1), reduced_data(idx_kmedoids == i, 2), reduced_data(idx_kmedoids == i, 3),...
        90, cluster_colors(i, :), 'filled');
    % scatter3(reduced_data(idx_kmedoids == i, 1), reduced_data(idx_kmedoids == i, 2), reduced_data(idx_kmedoids == i, 3),...
    % 50, cluster_colors(i, :), 'filled');
end
view(-45,45);
hold off;
title('K-Medoids Clustering (PCA-reduced Data)');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
grid on;
% Axis
ax=gca;
box off
axis square;
% Source Data
SourceDataTable = array2table([idx_kmedoids,reduced_data],'VariableNames',{'cluster idx','PC1','PC2','PC3'});
writetable(SourceDataTable,fullfile(output_dir, 'SourceData_KmedoidsClustering_PCA.csv'));
% Save the figure
[annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(output_dir, sprintf(['KmedoidsClustering_PC_visualization_clusters_clusters_' num2str(num_clusters) '_features_' num2str(length(feature_names)) '.pdf'])));
close(fig);

%% 📌 Display Results
fprintf('K-Means Silhouette Score: %.3f\n', silhouette_kmeans);
% fprintf('GMM Silhouette Score: %.3f\n', silhouette_gmm);
fprintf('K-Medoids Silhouette Score: %.3f\n', silhouette_kmedoids);

% %% 📌 Define Mixed-Strategy Mice
% mixed_threshold = 0.45;  % Define range for mixed classification
% mixed_mice_idx = (gm_probabilities(:,1) > mixed_threshold & gm_probabilities(:,1) < (1 - mixed_threshold));
% 
% % Count Mixed Mice
% num_mixed_mice = sum(mixed_mice_idx);
% fprintf('Number of Mixed-Strategy Mice (GMM): %d\n', num_mixed_mice);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract K-Means Cluster Labels
% Generate cluster indices dynamically
cluster_idx = cell(1, num_clusters);
for k = 1:num_clusters
    cluster_idx{k} = (idx_kmeans == k); % Logical indexing for each cluster
end

%% Extract Feature Data by Cluster
% Initialize an empty structure
cluster_features = struct();

% Loop through each feature name
for i = 1:length(feature_names)
    feature_name = feature_names{i}; % Get feature name as string
    
    % Extract feature data for each cluster and store in a cell array
    cluster_data = cell(1, num_clusters);
    for k = 1:num_clusters
        cluster_data{k} = features(cluster_idx{k}, i);
    end
    
    % Assign to structure
    cluster_features.(feature_name) = cluster_data;
end


%% Perform Statistical Comparisons (Permutation Tests)
clear pval Factual

% ANOVA for all features across clusters
for i = 1:size(features,2)
    [pval(i), Factual(i), ~] = randanova1(features(:,i), idx_kmeans, 1000);
end

% Define comparison pairs (Cluster 1 vs 2, 1 vs 3, 2 vs 3)
comparison_pairs = [1, 2; 1, 3; 2, 3];

% Initialize storage for results
num_comparisons = size(comparison_pairs, 1);
p_values = struct();
diff_obs_values = struct();
feature_names = fieldnames(cluster_features);

% Set permutation parameters
n_perm = 100;
method = 'median'; % Can be 'mean' or 'median'

% Perform permutation tests
for f = 1:length(feature_names)
    feature = feature_names{f};
    p_values.(feature) = zeros(1, num_comparisons);
    diff_obs_values.(feature) = zeros(1, num_comparisons);

    % Calculate real swap ratio
    % Initialize storage for within- and between-animal variability
    within_variability_values = zeros(1, 3);
    between_variability_values = zeros(1, 3);

    % Loop through all clusters (1, 2, and 3)
    for cluster = 1:num_clusters
        % Get unique animal IDs in the current cluster
        unique_animals = unique(id(cluster_idx{cluster}));

        % Compute within-animal variability (using min-max range)
        animal_variability = arrayfun(@(a) ...
            std(cluster_features.(feature){cluster}(strcmp(id(cluster_idx{cluster}), unique_animals{a}))), ...
            1:length(unique_animals));

        % Compute between-animal variability (range of per-animal medians)
        animal_means = arrayfun(@(a) ...
            mean(cluster_features.(feature){cluster}(strcmp(id(cluster_idx{cluster}), unique_animals{a}))), ...
            1:length(unique_animals));

        % Store median values for each cluster
        within_variability_values(cluster) = mean(animal_variability);
        between_variability_values(cluster) = std(animal_means);
    end

    % Estimate swap ratio
    % estimated_swap_ratio = between_variability / (within_variability + between_variability);
    estimated_swap_ratio = between_variability_values ./ (within_variability_values + between_variability_values);
    % estimated_swap_ratio = .5;
    disp(['Estimated Realistic Swap Ratio: ', num2str(estimated_swap_ratio)]);

    for c = 1:num_comparisons
        cluster_a = comparison_pairs(c, 1);
        cluster_b = comparison_pairs(c, 2);

        % Dynamically select the correct indices for each cluster
        idx_a = cluster_idx{cluster_a};
        idx_b = cluster_idx{cluster_b};

        mode='trial';

        [p_values.(feature)(c), diff_obs_values.(feature)(c)] = mixed_permutation_test_unpaired( ...
            cluster_features.(feature){cluster_a}, cluster_features.(feature){cluster_b}, ...
            id(idx_a), id(idx_b), n_perm, mean(estimated_swap_ratio([cluster_a,cluster_b])), method, mode, true);

        [h_values_KS.(feature)(c), p_values_KS.(feature)(c)] = kstest2( ...
            cluster_features.(feature){cluster_a}, cluster_features.(feature){cluster_b});
    end
end

%% Boxplots to Compare Clusters
fig = figure('Visible','on');
fig.Position = [100, 100, 1000, 450]; % [x, y, width, height]
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for nb_idx = 1:size(features,2)
    nexttile;
    feature = feature_names{nb_idx};
    plot_data = [];
    grouping_data = [];
    for cluster = 1:num_clusters
        plot_data = [plot_data;cluster_features.(feature){cluster}];
        grouping_data = [grouping_data;cluster*ones(size(cluster_features.(feature){cluster}))];
    end
    h = boxplot(plot_data, grouping_data, 'Symbol', 'o', 'OutlierSize', 0.01, 'Whisker', 1.5, 'colors', cluster_colors, 'boxstyle', 'outline');
    % Find and modify the properties of different boxplot elements
    set(h(1,:),'LineWidth',1.5,'LineStyle','-'); % Upper Whisker
    set(h(2,:),'LineWidth',1.5,'LineStyle','-'); % Lower Whisker
    set(h(3,:),'LineWidth',1.5,'LineStyle','-'); % Upper Adjacent Value
    set(h(4,:),'LineWidth',1.5,'LineStyle','-'); % Lower Adjacent Value
    set(h(5,:),'LineWidth',1.5,'LineStyle','-'); % Box
    set(h(6,:),'LineWidth',2.5,'LineStyle','-'); % Median

    
    hold on;

    % Overlay individual data points with colors matching cluster groups
    for cl = unique(grouping_data)'
        idx = (grouping_data == cl); % Find indices for each cluster
        sw{cl} = swarmchart(grouping_data(idx), plot_data(idx), 10, ...
            'filled', 'MarkerFaceColor', cluster_colors(cl, :), ...
            'MarkerFaceAlpha', 0.7, 'XJitter', 'density', 'XJitterWidth', 0.4);
    end
    
    % Set plot properties
    title(feature,'Interpreter','none','FontSize',14);

    % Cluster naming
    for i=1:num_clusters
        mean_hits(i)=mean(features(idx_kmeans==i,strcmp(feature_names,'correct_hits')));
        mean_rejections(i)=mean(features(idx_kmeans==i,strcmp(feature_names,'correct_rejections')));
    end
    [~,idx_hits]=min(mean_hits);
    cluster_names{idx_hits} = 'cautious';
    [~,idx_rej]=min(mean_rejections);
    cluster_names{idx_rej} = 'impulsive';
    cluster_names{~ismember([1,2,3],[idx_hits,idx_rej])} = 'balanced';

    xticklabels(cluster_names);
    % Axis
    ax=gca;
    ax.FontSize = 14;
    ax.LineWidth = 2;
    if contains(feature,'switch')
        ax.YLim(2)=300;
    end
    ax.YLim(1)=0;
    box off

    %% Display statistical results in the figure
    stats_text = {
        ['Permutation Results:']
        ['p1to2 = ' num2str(p_values.(feature)(1))]
        ['p1to3 = ' num2str(p_values.(feature)(2))]
        ['p2to3 = ' num2str(p_values.(feature)(3))]
        };

    % Position the text in the upper-right of the figure
    xlims = xlim;
    ylims = ylim;
    text(xlims(2) * 0.8, ylims(2) * 0.9, stats_text, 'FontSize', 8, 'BackgroundColor', 'w', 'Interpreter', 'none');

    % Source Data
    SourceDataTable = array2table([plot_data, grouping_data],'VariableNames',{feature ,'cluster idx',});
    writetable(SourceDataTable,fullfile(output_dir, ['SourceData_Boxplot_' feature '.csv']));
end

% title
sg = sgtitle('Behavioral Differences Between K-Means Clusters');
sg.Interpreter = 'none';

% Ensure font sizes and aspect ratios remain fixed
set(findall(fig, '-property', 'FontSize'), 'FontSize', 12); % Set all font sizes
set(gca, 'LooseInset', get(gca, 'TightInset')); % Remove extra padding

% Save the figure
[annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(output_dir, sprintf(['NotBoxPlot_BehavioralDifferences_between_KmeanClusters_clusters_' num2str(num_clusters) '_features_' num2str(length(feature_names)) '.pdf'])),'ContentType', 'vector', 'Resolution', 300); % Use vector graphics for no scaling issues
close(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Cluster Assignments to summary_data.lickport.D1_End
% List of variables to check and delete
fields = fieldnames(summary_data.lickport.D1_End);
vars_to_remove = fields(contains(fields,{'kmeans_cluster', 'gmm_cluster', 'gmm_probability', 'kmedoids_cluster',}));

% Check if the field exists and remove it
for i = 1:length(vars_to_remove)
    summary_data.lickport.D1_End = removevars(summary_data.lickport.D1_End, vars_to_remove{i});
    fprintf('Removed field: %s\n', vars_to_remove{i});
end
% Assign Kmeans Clusters, GMM Clusters, GMM Probabilities (posterior probability of belonging to Cluster 1)
clusterResultsTable = table(id,cohort,idx_kmeans,idx_kmedoids, ...
    'VariableNames',{'Mouse_RFID','cohort','kmeans_cluster','kmedoids_cluster'});
% Join tables
summary_data.lickport.D1_End = outerjoin(summary_data.lickport.D1_End, clusterResultsTable, ...
    'Keys', {'Mouse_RFID', 'cohort'}, ...
    'MergeKeys', true, ...   % Keep matching keys together
    'LeftVariables', summary_data.lickport.D1_End.Properties.VariableNames, ... % Keep all original columns
    'RightVariables', setdiff(clusterResultsTable.Properties.VariableNames, {'Mouse_RFID', 'cohort'}) ... % Only add non-key columns
    );

%% Save the updated summary_data structure (optional)
save(fullfile(main_dir, 'data', 'processed', 'cross_cohort_files', 'summary_data.mat'), 'summary_data');

%% Display Results
fprintf('K-Means and GMM clusters successfully added to summary_data.lickport.D1_End\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis of association of clusering and social hierarchy and chasing
% Clear previous data
clear myTable

% Merge relevant datasets based on Mouse RFID and cohort information
myTable = join(summary_data.lickport.D1_End, summary_data.tube.D1_End, 'Keys', {'Mouse_RFID','cohort'});

% Filter out rows where k-means clustering data is NaN
myTable = myTable(~isnan(myTable.kmeans_cluster), :);

% Define "social" metrics for analysis
myMetrics = {'DSz_Competition', ...
    'Fraction_Wins_Competition', 'Fraction_Losses_Competition', ...
    'N_Events_Chasing', ...
    'Fraction_Wins_Chasing', 'Fraction_Losses_Chasing'};

%% Improved Boxplots with NotBoxPlot
clear p_values diff_obs_values

fig = figure('Visible','on');
fig.Position = [100, 100, 1000, 450]; % [x, y, width, height]
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

numMetrics = length(myMetrics);
numClusters = length(unique(myTable.kmeans_cluster));

% Define comparison pairs (Cluster 1 vs 2, 1 vs 3, 2 vs 3)
comparison_pairs = [1, 2; 1, 3; 2, 3];
% Initialize storage for results
num_comparisons = size(comparison_pairs, 1);
p_values = struct();
diff_obs_values = struct();
p_values = struct;
% Set permutation parameters
n_perm = 1000;
method = 'median'; % Can be 'mean' or 'median'

% Loop through each metric and plot using notBoxPlot
for ix = 1:numMetrics
    nexttile;
    plot_data = myTable.(myMetrics{ix});
    grouping_data = myTable.kmeans_cluster;
    current_metric = myMetrics{ix};

    % nb{ix} = notBoxPlot(metricData, clusterGroups, 'jitter', 0.5, 'markMedian', true);
    % Boxplot visualization
    h = boxplot(plot_data, grouping_data, 'Symbol', 'o', 'OutlierSize', 0.01, 'Whisker', 1.5, 'colors', cluster_colors, 'boxstyle', 'outline');
    % Find and modify the properties of different boxplot elements
    set(h(1,:),'LineWidth',1.5,'LineStyle','-'); % Upper Whisker
    set(h(2,:),'LineWidth',1.5,'LineStyle','-'); % Lower Whisker
    set(h(3,:),'LineWidth',1.5,'LineStyle','-'); % Upper Adjacent Value
    set(h(4,:),'LineWidth',1.5,'LineStyle','-'); % Lower Adjacent Value
    set(h(5,:),'LineWidth',1.5,'LineStyle','-'); % Box
    set(h(6,:),'LineWidth',2.5,'LineStyle','-'); % Median
    hold on;

    % Overlay individual data points with colors matching cluster groups
    for cl = unique(grouping_data)'
        idx = (grouping_data == cl); % Find indices for each cluster
        sw{cl} = swarmchart(grouping_data(idx), plot_data(idx), 10, ...
            'filled', 'MarkerFaceColor', cluster_colors(cl, :), ...
            'MarkerFaceAlpha', 0.7, 'XJitter', 'density', 'XJitterWidth', 0.4);
    end
    hold on;

    % Set plot properties
    title(myMetrics{ix}, 'Interpreter', 'none', 'FontSize', 14);
    xticklabels(cluster_names);

    % Customize plot aesthetics
    ax = gca;
    ax.FontSize = 14;
    ax.LineWidth = 2;
    box off

    %% Statistical Analysis: Permutation Test
    % Calculate real swap ratio
    % Initialize storage for within- and between-animal variability
    within_variability_values = zeros(1, 3);
    between_variability_values = zeros(1, 3);

    % Loop through all clusters (1, 2, and 3)
    for cluster = 1:num_clusters
        % Get unique animal IDs in the current cluster
        unique_animals = unique(id(cluster_idx{cluster}));

        % Compute within-animal variability (using min-max range)
        animal_variability = arrayfun(@(a) ...
            std(myTable.(current_metric)(strcmp(id(cluster_idx{cluster}), unique_animals{a}))), ...
            1:length(unique_animals));

        % Compute between-animal variability (range of per-animal medians)
        animal_means = arrayfun(@(a) ...
            mean(myTable.(current_metric)(strcmp(id(cluster_idx{cluster}), unique_animals{a}))), ...
            1:length(unique_animals));

        % Store median values for each cluster
        within_variability_values(cluster) = mean(animal_variability);
        between_variability_values(cluster) = std(animal_means);
    end

    % Aggregate across clusters
    % within_variability = mean(within_variability_values); % Use median to avoid outlier effects
    % between_variability = mean(between_variability_values); % Same for between-animal variability

    % Estimate swap ratio
    % estimated_swap_ratio = between_variability / (within_variability + between_variability);
    estimated_swap_ratio = between_variability_values ./ (within_variability_values + between_variability_values);
    % estimated_swap_ratio = .5;
    disp(['Estimated Realistic Swap Ratio: ', num2str(estimated_swap_ratio)]);

    for c = 1:num_comparisons
        cluster_a = comparison_pairs(c, 1);
        cluster_b = comparison_pairs(c, 2);

        % Dynamically select the correct indices for each cluster
        idx_a = cluster_idx{cluster_a};
        idx_b = cluster_idx{cluster_b};

        [p_values.(current_metric)(c), diff_obs_values.(current_metric)(c)] = mixed_permutation_test_unpaired( ...
            myTable.(current_metric)(idx_a), myTable.(current_metric)(idx_b), ...
            id(idx_a), id(idx_b), n_perm, mean(estimated_swap_ratio([cluster_a,cluster_b])), method, 'trial', true);


        [h_values_KS.(current_metric)(c), p_values_KS.(current_metric)(c)] = kstest2( ...
            myTable.(current_metric)(idx_a), myTable.(current_metric)(idx_b));
    end

    % Display Statistical Results in the Figure
    hold on;
    stats_text = {
        'Permutation Results:', ...
        ['p1to2 = ', num2str(p_values.(myMetrics{ix})(1))], ...
        ['p1to3 = ', num2str(p_values.(myMetrics{ix})(2))], ...
        ['p2to3 = ', num2str(p_values.(myMetrics{ix})(3))]
        };

    % Position text in the upper-right of the subplot
    xlims = xlim;
    ylims = ylim;
    text(xlims(2) * 0.8, ylims(2) * 0.9, stats_text, 'FontSize', 8, ...
        'BackgroundColor', 'w', 'Interpreter', 'none');

    % Source Data
    SourceDataTable = array2table([plot_data, grouping_data],'VariableNames',{current_metric ,'cluster idx',});
    writetable(SourceDataTable,fullfile(output_dir, ['SourceData_Boxplot_' current_metric '.csv']));
end

%% Add Figure Title and Save Output
sgtitle('Social Traits Differences Between K-Means Clusters', 'Interpreter', 'none');

% Ensure font sizes and aspect ratios remain fixed
set(findall(fig, '-property', 'FontSize'), 'FontSize', 12); % Set all font sizes
set(gca, 'LooseInset', get(gca, 'TightInset')); % Remove extra padding
hold off

% Save the figure
[annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
output_filepath = fullfile(output_dir, sprintf('NotBoxPlot_SocialTraitsDifferences_between_KmeanClusters_clusters_%d_features_%d.pdf', numClusters, length(feature_names)));
exportgraphics(fig, output_filepath);
close(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Repetitions Analysis
clear myTable
current_data = summary_data.lickport.D1_End;

% Get unique animals and their repetitions
unique_animals = unique(current_data.Mouse_RFID);
repetitions_index = cell(size(unique_animals));

for an = 1:numel(unique_animals)
    animal_rows = find(strcmp(current_data.Mouse_RFID, unique_animals{an}));
    repetitions_index{an} = animal_rows;
end

%% Analyze Stability for Repetitions
ii = 1; % First repetition
jj = 2; % Second repetition

all_data = nan(length(unique_animals),max(cellfun(@numel,repetitions_index)));
counter = 1;

for an = 1:numel(unique_animals)
    for rep = 1:length(repetitions_index{an})
        all_data(an,rep)=current_data.kmeans_cluster(repetitions_index{an}(rep));
    end
    % if numel(repetitions_index{an}) >= jj
    %     xData = cat(1, xData, current_data.(metric)(repetitions_index{an}(ii)));
    %     yData = cat(1, yData, current_data.(metric)(repetitions_index{an}(jj)));
    %     animal_ID{counter,1} = unique_animals{an};
    selection_idx(an) = repetitions_index{an}(ii);
    %     counter = counter + 1;
    % end
end

% Apply function to move non-NaN values to the beginning of each row,
% while keeping the order of appearance, and pushing NaNs to the end.
sortedMatrix = sortNonNaNFirst(all_data);

% Data
valid_idx = ~isnan(sortedMatrix(:,ii)) & ~isnan(sortedMatrix(:,jj));
xData = sortedMatrix(valid_idx,ii);
yData = sortedMatrix(valid_idx,jj);
animal_ID = unique_animals(valid_idx);
selection_idx = selection_idx(valid_idx);
xData_named = cluster_names(xData);
yData_named = cluster_names(yData);

myTable = table(animal_ID, categorical(xData_named'), categorical(yData_named'), 'VariableNames', {'Mouse_RFID', 'kmeans_cluster_R1', 'kmeans_cluster_R2'});
writetable(myTable,fullfile(output_dir, ['SourceData_RepetitionAnalysis.csv']));

%% Sankey Diagram

% Use table for sankey flow chart
myData = myTable(:,2:end);

% Customizable options
% Colormap: can be the name of matlab colormaps or a matrix of (N x 3).
%   Important: N must be the max number of categories in a layer 
%   multiplied by the number of layers. 
%   In the example of Option 1, N should be 4 * 4 = 16;
%   In the example of Option 2, N should be 5 * 3 = 15;
options.color_map = 'jet';      
options.flow_transparency = 0.2;   % opacity of the flow paths
options.bar_width = 40;            % width of the category blocks
options.show_perc = false;          % show percentage over the blocks
options.text_color = [0 0 0];      % text color for the percentages
options.show_layer_labels = true;  % show layer names under the chart
options.show_cat_labels = true;   % show categories over the blocks.
options.show_legend = false;        % show legend with the category names. 
                                   % if the data is not a table, then the
                                   % categories are labeled as catX-layerY
plotSankeyFlowChart(myData,options);
% Save the figure
[annot, srcInfo] = docDataSrc(gcf,fullfile(output_dir),mfilename('fullpath'),logical(1));
exportgraphics(gcf, fullfile(output_dir, sprintf(['SankeyDiagram_StabilityOverRounds_Clusters' num2str(num_clusters) '.pdf'])));
close(gcf);

%% Contingency matrix
fig=figure;
% Define unique values in both datasets
unique_x = unique(xData);
unique_y = unique(yData);

% Create contingency matrix (initialize with zeros)
contingency_matrix = zeros(length(unique_x), length(unique_y));

% Count occurrences of each (x, y) pair
for i = 1:length(xData)
    x_idx = find(unique_x == xData(i));  % Find index for x
    y_idx = find(unique_y == yData(i));  % Find index for y
    contingency_matrix(y_idx, x_idx) = contingency_matrix(y_idx, x_idx) + 1;
end

% plot
imagesc(contingency_matrix);    
crameri('vik');
colorbar;
xticks(1:3);yticks(1:3);
xticklabels(cluster_names);
yticklabels(cluster_names);
xlabel('round 1');
ylabel('round 2');
 
% Overlay text with the actual values in the matrix
[num_rows, num_cols] = size(contingency_matrix);
for row = 1:num_rows
    for col = 1:num_cols
        text(col, row, num2str(contingency_matrix(row, col)), ...
            'FontSize', 12, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'Color', [1,1,1]);
    end
end

[annot, srcInfo] = docDataSrc(fig,fullfile(output_dir),mfilename('fullpath'),logical(1));
exportgraphics(fig, fullfile(output_dir, sprintf(['ContingencyMatrix_StabilityOverRounds_Clusters' num2str(num_clusters) '.pdf'])));
close(fig);

% Check WT/OXtKO in the clusters
myData = summary_data.lickport.D1_End;
valid_idx = ~isnan(myData.kmeans_cluster);
genotype_clean = myData.genotype_corentin(valid_idx);
cluster_clean = myData.kmeans_cluster(valid_idx);
[tbl, chi2stat, p, labels] = crosstab(genotype_clean, cluster_clean);
row_labels = labels(1:2,1);  % WT, OxtKO
col_labels = labels(1:3,2);  % Cluster 1, 2, 3
disp(array2table(tbl, 'RowNames', row_labels, 'VariableNames', strcat("Cluster", col_labels)));
fprintf('Chi-square test p-value: %.4f\n', p);


%%
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