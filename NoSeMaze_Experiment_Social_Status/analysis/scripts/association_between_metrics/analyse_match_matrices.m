%% plot_match_matrices_tube_metrics.m
% Author: Jonathan Reinwald
% Date: January 2025
% Purpose: Process match matrices, compute David's Scores, and visualize
% metrics for tube competition and chasing across cohorts and day ranges.

%% Pre-Clearing
clear all; clc; close all;

%% Set Directories
main_dir = 'myRootPath/NoSeMaze_Experiment'; % --> NOTE: Replace with own root directory
processed_dir = fullfile(main_dir, 'data', 'processed');
results_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort', 'associations_between_metrics', 'tube', 'match_matrices');
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');

% Add helper paths
addpath(genpath(fullfile(main_dir, 'src')));

% Ensure the results directory exists
if ~isfolder(results_dir)
    mkdir(results_dir);
end

%% Cohort Selection
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));
cohortsTbl = cohortsTbl(cohortsTbl.use_tube == 1, :);

%% Parameters
day_ranges = {'D1_21'};%,'D1_14'};%
day_ranges_idx = {[1, 21]};%,[1, 14]};%
normalizations = {'maximum'}; %{'none', 'maximum', 'InteractionInDyad', 'events'};
chasing_sortingType = 'FractionWins'; % 'FractionWins'DSz_Chasing

%% Processing and Visualization
% Loop through each normalization method
for norm_idx = 1:length(normalizations)
    normalization = normalizations{norm_idx};

    % Loop through each day range
    for day_idx = 1:length(day_ranges)
        day_range = day_ranges{day_idx};
        day_range_idx = day_ranges_idx{day_idx};

        % Reinitialize match matrix container
        myMatchMatrix = [];

        % Loop through each cohort
        for cohort_idx = 1:height(cohortsTbl)
            cohort = cohortsTbl.cohort{cohort_idx};

            % Load match matrix data
            cohort_data_dir = fullfile(processed_dir, cohort, 'tube', 'full_hierarchy_files');
            full_hierarchy_file = fullfile(cohort_data_dir, ['full_hierarchy_' cohort '.mat']);

            if ~isfile(full_hierarchy_file)
                warning('File not found for cohort: %s', cohort);
                continue;
            end

            load(full_hierarchy_file, 'full_hierarchy');

            % Compute match matrices
            current_match_matrix_Competition = sum(cat(3, full_hierarchy(day_range_idx(1):day_range_idx(2)).match_matrix), 3);
            current_match_matrix_Chasing = sum(cat(3, full_hierarchy(day_range_idx(1):day_range_idx(2)).match_matrix_chasing), 3);

            % Compute David's Scores
            DS_info = compute_DS_from_match_matrix(current_match_matrix_Competition);
            myData.Competition.DSz{cohort_idx} = zscore(DS_info.DS)';
            DS_info = compute_DS_from_match_matrix(current_match_matrix_Chasing);
            myData.Chasing.DSz{cohort_idx} = zscore(DS_info.DS)';
            % myData.Competition.FractionWins{cohort_idx} = sum(current_match_matrix_Competition, 2)./sum(current_match_matrix_Competition(:));
            % myData.Competition.FractionLosses{cohort_idx} = sum(current_match_matrix_Competition, 1)./sum(current_match_matrix_Competition(:));
            myData.Chasing.FractionWins{cohort_idx} = sum(current_match_matrix_Chasing, 2)./sum(current_match_matrix_Chasing(:));
            myData.Chasing.FractionLosses{cohort_idx} = sum(current_match_matrix_Chasing, 1)./sum(current_match_matrix_Chasing(:));

            % Sorting indices
            [~, sorting_idx_DSzChas] = sort(myData.Chasing.DSz{cohort_idx}, 'descend');
            [~, sorting_idx_DSzComp] = sort(myData.Competition.DSz{cohort_idx}, 'descend');
            [~, sorting_idx_FractionChas] = sort(myData.Chasing.FractionWins{cohort_idx}, 'descend');

            % Apply normalization and store sorted matrices
            if strcmp(chasing_sortingType,'FractionWins')
                myMatchMatrix = apply_normalization(myMatchMatrix, current_match_matrix_Competition, ...
                    current_match_matrix_Chasing, sorting_idx_DSzComp, sorting_idx_FractionChas, cohort_idx, normalization);
            elseif strcmp(chasing_sortingType,'DSz_Chasing')
                myMatchMatrix = apply_normalization(myMatchMatrix, current_match_matrix_Competition, ...
                    current_match_matrix_Chasing, sorting_idx_DSzComp, sorting_idx_DSzChas, cohort_idx, normalization);
            end
        end


        %% Quantification & Visualization of Quadrant Proportions

        % 1) Determine the number of cohorts dynamically from your data structure.
        %    We assume 'Chasing_sortedByComp' and 'Competition_sortedByComp' both
        %    have the same number of cohorts.
        numCohorts = numel(myMatchMatrix.Chasing_sortedByComp);

        % Preallocate vectors for quadrant proportions (chasing)
        chasingbycompetition_topleft = nan(numCohorts,1);
        chasingbycompetition_topright = nan(numCohorts,1);
        chasingbycompetition_bottomleft = nan(numCohorts,1);
        chasingbycompetition_bottomright = nan(numCohorts,1);

        chasingbychasing_topleft = nan(numCohorts,1);
        chasingbychasing_topright = nan(numCohorts,1);
        chasingbychasing_bottomleft = nan(numCohorts,1);
        chasingbychasing_bottomright = nan(numCohorts,1);

        % Preallocate vectors for quadrant proportions (competition)
        competition_topleft = nan(numCohorts,1);
        competition_topright = nan(numCohorts,1);
        competition_bottomleft = nan(numCohorts,1);
        competition_bottomright = nan(numCohorts,1);

        top_down_chasing = nan(numCohorts,1);
        bottom_up_chasing = nan(numCohorts,1);
        down_vs_up_chasing = nan(numCohorts,1);
        top_down_chasing_weighted = nan(numCohorts,1);
        bottom_up_chasing_weighted = nan(numCohorts,1);
        down_vs_up_chasing_weighted = nan(numCohorts,1);

        chasingbycompetition_topleft_all = nan(numCohorts,25);
        chasingbycompetition_topright_all = nan(numCohorts,25);
        chasingbycompetition_bottomleft_all = nan(numCohorts,25);
        chasingbycompetition_bottomright_all = nan(numCohorts,25);

        chasingbychasing_topleft_all = nan(numCohorts,25);
        chasingbychasing_topright_all = nan(numCohorts,25);
        chasingbychasing_bottomleft_all = nan(numCohorts,25);
        chasingbychasing_bottomright_all = nan(numCohorts,25);

        competition_topleft_all = nan(numCohorts,25);
        competition_topright_all = nan(numCohorts,25);
        competition_bottomleft_all = nan(numCohorts,25);
        competition_bottomright_all = nan(numCohorts,25);

        %% 2) Calculate quadrant proportions for each cohort
        for ix = 1:numCohorts
            % Extract the current chasing matrix (sorted by competition DS)
            chasingbycompetitionMat = myMatchMatrix.Chasing_sortedByComp{ix};

            % Quadrant sums: top-left (1:5,1:5), top-right (1:5,6:10),
            %                bottom-left (6:10,1:5), bottom-right (6:10,6:10)
            chasingbycompetitionMat(1:size(chasingbycompetitionMat,1)+1:end) = NaN;
            totalChasing = nansum(nansum(chasingbycompetitionMat(1:10,1:10)));

            %
            top_down_chasing(ix) = nansum(nansum(triu(chasingbycompetitionMat)))./nansum(nansum(chasingbycompetitionMat));
            bottom_up_chasing(ix) = nansum(nansum(tril(chasingbycompetitionMat)))./nansum(nansum(chasingbycompetitionMat));
            down_vs_up_chasing(ix) =  top_down_chasing(ix)./bottom_up_chasing(ix);

            % Weighting matrix
            % Top-Down
            n = 10;  % number of animals
            % Get matrices of row and column indices
            [rows, cols] = ndgrid(1:n, 1:n);
            % Rank difference (positive = top-down)
            rankDiff = cols - rows;
            % Assign weight: from 9 (for diff = 1) down to 1 (for diff = 9)
            weights_TopDown = (n - 1) - rankDiff + 1;  % e.g., for rankDiff = 1 → weight = 9
            weights_TopDown(rankDiff <= 0) = 0;       % only consider top-down
            % Bottom-Up
            [rows, cols] = ndgrid(1:n, 1:n);
            rankDiff = rows - cols;  % i > j = bottom-up
            weights_BottomUp = rankDiff;      % larger diff = higher weight
            weights_BottomUp(rankDiff <= 0) = 0;  % keep only bottom-up interactions
            weights = weights_BottomUp./sum(sum(weights_BottomUp))+weights_TopDown./sum(sum(weights_TopDown));
            % % n = 10;
            % % [rows, cols] = ndgrid(1:n, 1:n);
            % % % Top-down weights: i < j → higher weight for closer ranks
            % % rankDiff_TD = cols - rows;
            % % weights_TopDown = zeros(n);
            % % weights_TopDown(rankDiff_TD > 0) = 1 ./ rankDiff_TD(rankDiff_TD > 0);
            % % % Bottom-up weights: i > j → higher weight for *larger* gaps
            % % rankDiff_BU = rows - cols;
            % % weights_BottomUp = zeros(n);
            % % weights_BottomUp(rankDiff_BU > 0) = 1 ./ (n - rankDiff_BU(rankDiff_BU > 0));
            % % % -------- COMBINED --------
            % % weights = weights_TopDown + weights_BottomUp;
            %
            top_down_chasing_weighted(ix) = nansum(nansum(triu(chasingbycompetitionMat).*weights))./nansum(nansum(chasingbycompetitionMat.*weights));
            bottom_up_chasing_weighted(ix) = nansum(nansum(tril(chasingbycompetitionMat).*weights))./nansum(nansum(chasingbycompetitionMat.*weights));
            down_vs_up_chasing_weighted(ix) =  top_down_chasing_weighted(ix)./bottom_up_chasing_weighted(ix);

            % chasingbycompetition_topleft(ix)     = nansum(nansum(chasingbycompetitionMat(1:5,  1:5)))  / totalChasing;
            % chasingbycompetition_topright(ix)    = nansum(nansum(chasingbycompetitionMat(1:5,  6:10))) / totalChasing;
            % chasingbycompetition_bottomleft(ix)  = nansum(nansum(chasingbycompetitionMat(6:10, 1:5)))  / totalChasing;
            % chasingbycompetition_bottomright(ix) = nansum(nansum(chasingbycompetitionMat(6:10, 6:10))) / totalChasing;

            % Get sums and non-NaN counts for each quadrant
            sum_topleft     = nansum(chasingbycompetitionMat(1:5,  1:5), 'all');
            sum_topright    = nansum(chasingbycompetitionMat(1:5,  6:10), 'all');
            sum_bottomleft  = nansum(chasingbycompetitionMat(6:10, 1:5), 'all');
            sum_bottomright = nansum(chasingbycompetitionMat(6:10, 6:10), 'all');

            count_topleft     = sum(~isnan(chasingbycompetitionMat(1:5,  1:5)), 'all');
            count_topright    = sum(~isnan(chasingbycompetitionMat(1:5,  6:10)), 'all');
            count_bottomleft  = sum(~isnan(chasingbycompetitionMat(6:10, 1:5)), 'all');
            count_bottomright = sum(~isnan(chasingbycompetitionMat(6:10, 6:10)), 'all');

            % Compute means (per-entry chasing)
            mean_topleft     = sum_topleft     / count_topleft;
            mean_topright    = sum_topright    / count_topright;
            mean_bottomleft  = sum_bottomleft  / count_bottomleft;
            mean_bottomright = sum_bottomright / count_bottomright;

            % Normalize so the four mean values sum to 1
            total_mean = mean_topleft + mean_topright + mean_bottomleft + mean_bottomright;

            chasingbycompetition_topleft(ix)     = mean_topleft     / total_mean;
            chasingbycompetition_topright(ix)    = mean_topright    / total_mean;
            chasingbycompetition_bottomleft(ix)  = mean_bottomleft  / total_mean;
            chasingbycompetition_bottomright(ix) = mean_bottomright / total_mean;

            vec = chasingbycompetitionMat(1:5,1:5);
            vec = vec(~isnan(vec));
            chasingbycompetition_topleft_all(ix,1:length(vec)) = vec;
            vec = chasingbycompetitionMat(1:5,6:10);
            vec = vec(~isnan(vec));
            chasingbycompetition_topright_all(ix,1:length(vec)) = vec;
            vec = chasingbycompetitionMat(6:10,1:5);
            vec = vec(~isnan(vec));
            chasingbycompetition_bottomleft_all(ix,1:length(vec))  = vec;
            vec = chasingbycompetitionMat(6:10,6:10);
            vec = vec(~isnan(vec));
            chasingbycompetition_bottomright_all(ix,1:length(vec)) = vec;

            % Extract the current competition matrix (sorted by competition DS)
            compMat = myMatchMatrix.Competition_sortedByComp{ix};%Competition_sortedByComp
            compMat(1:size(compMat,1)+1:end) = NaN;
            % Get sums and non-NaN counts for each quadrant
            sum_topleft     = nansum(compMat(1:5,  1:5), 'all');
            sum_topright    = nansum(compMat(1:5,  6:10), 'all');
            sum_bottomleft  = nansum(compMat(6:10, 1:5), 'all');
            sum_bottomright = nansum(compMat(6:10, 6:10), 'all');

            count_topleft     = sum(~isnan(compMat(1:5,  1:5)), 'all');
            count_topright    = sum(~isnan(compMat(1:5,  6:10)), 'all');
            count_bottomleft  = sum(~isnan(compMat(6:10, 1:5)), 'all');
            count_bottomright = sum(~isnan(compMat(6:10, 6:10)), 'all');

            % Compute means (per-entry chasing)
            mean_topleft     = sum_topleft     / count_topleft;
            mean_topright    = sum_topright    / count_topright;
            mean_bottomleft  = sum_bottomleft  / count_bottomleft;
            mean_bottomright = sum_bottomright / count_bottomright;

            % Normalize so the four mean values sum to 1
            total_mean = mean_topleft + mean_topright + mean_bottomleft + mean_bottomright;

            competition_topleft(ix)     = mean_topleft     / total_mean;
            competition_topright(ix)    = mean_topright    / total_mean;
            competition_bottomleft(ix)  = mean_bottomleft  / total_mean;
            competition_bottomright(ix) = mean_bottomright / total_mean;

            compMat(1:size(compMat,1)+1:end) = NaN;
            vec = compMat(1:5,1:5);
            vec = vec(~isnan(vec));
            competition_topleft_all(ix,1:length(vec)) = vec;
            vec = compMat(1:5,6:10);
            vec = vec(~isnan(vec));
            competition_topright_all(ix,1:length(vec)) = vec;
            vec = compMat(6:10,1:5);
            vec = vec(~isnan(vec));
            competition_bottomleft_all(ix,1:length(vec))  = vec;
            vec = compMat(6:10,6:10);
            vec = vec(~isnan(vec));
            competition_bottomright_all(ix,1:length(vec)) = vec;

            %% Extract the current chasing matrix (sorted by competition DS)
            chasingbychasingMat = myMatchMatrix.Chasing_sortedByChas{ix};
            chasingbychasingMat(1:size(chasingbychasingMat,1)+1:end) = NaN;

            % Get sums and non-NaN counts for each quadrant
            sum_topleft     = nansum(chasingbychasingMat(1:5,  1:5), 'all');
            sum_topright    = nansum(chasingbychasingMat(1:5,  6:10), 'all');
            sum_bottomleft  = nansum(chasingbychasingMat(6:10, 1:5), 'all');
            sum_bottomright = nansum(chasingbychasingMat(6:10, 6:10), 'all');

            count_topleft     = sum(~isnan(chasingbychasingMat(1:5,  1:5)), 'all');
            count_topright    = sum(~isnan(chasingbychasingMat(1:5,  6:10)), 'all');
            count_bottomleft  = sum(~isnan(chasingbychasingMat(6:10, 1:5)), 'all');
            count_bottomright = sum(~isnan(chasingbychasingMat(6:10, 6:10)), 'all');

            % Compute means (per-entry chasing)
            mean_topleft     = sum_topleft     / count_topleft;
            mean_topright    = sum_topright    / count_topright;
            mean_bottomleft  = sum_bottomleft  / count_bottomleft;
            mean_bottomright = sum_bottomright / count_bottomright;

            % Normalize so the four mean values sum to 1
            total_mean = mean_topleft + mean_topright + mean_bottomleft + mean_bottomright;

            chasingbychasing_topleft(ix)     = mean_topleft     / total_mean;
            chasingbychasing_topright(ix)    = mean_topright    / total_mean;
            chasingbychasing_bottomleft(ix)  = mean_bottomleft  / total_mean;
            chasingbychasing_bottomright(ix) = mean_bottomright / total_mean;

            vec = chasingbychasingMat(1:5,1:5);
            vec = vec(~isnan(vec));
            chasingbychasing_topleft_all(ix,1:length(vec)) = vec;
            vec = chasingbycompetitionMat(1:5,6:10);
            vec = vec(~isnan(vec));
            chasingbychasing_topright_all(ix,1:length(vec)) = vec;
            vec = chasingbycompetitionMat(6:10,1:5);
            vec = vec(~isnan(vec));
            chasingbychasing_bottomleft_all(ix,1:length(vec))  = vec;
            vec = chasingbycompetitionMat(6:10,6:10);
            vec = vec(~isnan(vec));
            chasingbychasing_bottomright_all(ix,1:length(vec)) = vec;
        end

        %% WRITE TABLE FOR R
        % Initialize
        data_ChasByComp = []; data_ChasByChas = []; data_CompByComp = [];
        quadrant = [];
        group_id = [];

        % Define quadrant names and corresponding matrices
        quadrant_names = ["Q1", "Q2", "Q3", "Q4"];
        matrices_ChasByComp = {
            chasingbycompetition_topleft_all, ...
            chasingbycompetition_topright_all, ...
            chasingbycompetition_bottomleft_all, ...
            chasingbycompetition_bottomright_all
            };
        matrices_ChasByChas = {
            chasingbychasing_topleft_all, ...
            chasingbychasing_topright_all, ...
            chasingbychasing_bottomleft_all, ...
            chasingbychasing_bottomright_all
            };
        matrices_CompByComp = {
            competition_topleft_all, ...
            competition_topright_all, ...
            competition_bottomleft_all, ...
            competition_bottomright_all
            };

        % Loop over quadrants
        for q = 1:4
            M = matrices_ChasByComp{q};  % Current matrix
            for g = 1:size(M, 1)  % For each group (row)
                vals = M(g, :);
                valid_idx = ~isnan(vals);
                data_ChasByComp = [data_ChasByComp; vals(valid_idx)'];
                quadrant = [quadrant; repmat(quadrant_names(q), sum(valid_idx), 1)];
                group_id = [group_id; repmat(g, sum(valid_idx), 1)];
            end
            M = matrices_CompByComp{q};  % Current matrix
            for g = 1:size(M, 1)  % For each group (row)
                vals = M(g, :);
                valid_idx = ~isnan(vals);
                data_CompByComp = [data_CompByComp; vals(valid_idx)'];
            end
            M = matrices_ChasByChas{q};  % Current matrix
            for g = 1:size(M, 1)  % For each group (row)
                vals = M(g, :);
                valid_idx = ~isnan(vals);
                data_ChasByChas = [data_ChasByChas; vals(valid_idx)'];
            end
        end
        
        T = table(data_ChasByComp, data_CompByComp, data_ChasByChas, categorical(quadrant), categorical(group_id), ...
            'VariableNames', {'ChasByComp', 'CompByComp', 'ChasByChas', 'quadrant', 'group_id'});
        writetable(T, fullfile(results_dir,'interaction_data.csv'));

        %% 3) Plot quadrant proportions for Chasing
        fig = figure('Name','Chasing Quadrant Proportions','Color','w','Visible','on','Position', [100, 100, 1000, 1400]);  % [left, bottom, width, height]

        cluster_colors = [0 0 255/2; 0 0 255/2; 0 0 255/2; 0 0 255/2; 0 0 255/2] / 255;

        subplot(3,2,1)

        % Combine the four quadrant vectors into columns for notBoxPlot
        chasingbycompetitionMatrix = [chasingbycompetition_topleft, chasingbycompetition_topright, chasingbycompetition_bottomleft, chasingbycompetition_bottomright];

        nb=notBoxPlot(chasingbycompetitionMatrix);
        for nb_idx = 1:length(nb)
            nb(nb_idx).data.MarkerEdgeColor = 'none';
            nb(nb_idx).data.MarkerFaceColor = cluster_colors(nb_idx,:);
            nb(nb_idx).data.MarkerSize = 3;
            nb(nb_idx).mu.Color = cluster_colors(nb_idx,:);
            nb(nb_idx).sdPtch(1).FaceColor = cluster_colors(nb_idx,:);
            nb(nb_idx).sdPtch(1).FaceAlpha = 0.2;
            nb(nb_idx).semPtch(1).FaceColor = cluster_colors(nb_idx,:);
            nb(nb_idx).semPtch(1).FaceAlpha = 0.4;
            nb(nb_idx).sdPtch(1).EdgeColor = 'none';
            % nb(nb_idx).sdPtch(2).EdgeColor = 'none';
            nb(nb_idx).semPtch(1).EdgeColor = 'none';
            % nb(nb_idx).semPtch(2).EdgeColor = 'none';
        end

        % Improve axes
        ax = gca;
        ax.XTick = 1:4;
        xlabels = {'Top-Left','Top-Right','Bottom-Left','Bottom-Right'};
        ax.XTickLabel = xlabels;
        ylabel('Proportion of total chasing');
        title('Chasing (sorted by Competition DS)');
        ylim([0 1]);

        hold on; % To add permutation-test results on the figure

        %% 4) Plot quadrant proportions for Chasing (high-resolution data)
        subplot(3,2,2)

        % Combine the four quadrant vectors into columns for notBoxPlot
        chasingbycompetitionMatrix_all = [chasingbycompetition_topleft_all(:), chasingbycompetition_topright_all(:), chasingbycompetition_bottomleft_all(:), chasingbycompetition_bottomright_all(:)];
        % chasingbycompetitionMatrix_all = [chasingbycompetition_topleft_all(:).^(1/3), chasingbycompetition_topright_all(:).^(1/3), chasingbycompetition_bottomleft_all(:).^(1/3), chasingbycompetition_bottomright_all(:).^(1/3)];
        % % % t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        % % nexttile;
        cluster_colors = [0 0 255/2; 0 0 255/2; 0 0 255/2; 0 0 255/2; 0 0 255/2] / 255;
        plot_data = [];
        grouping_data = [];
        for q_idx = 1:size(chasingbycompetitionMatrix_all,2)
            plot_data = [plot_data;chasingbycompetitionMatrix_all(:,q_idx)];
            grouping_data = [grouping_data;q_idx*ones(size(chasingbycompetitionMatrix_all,1),1)];
        end
        h = boxplot(plot_data, grouping_data, 'Symbol', 'o', 'OutlierSize', 0.01, 'Whisker', 1, 'colors', cluster_colors, 'boxstyle', 'outline');
        % Find and modify the properties of different boxplot elements
        set(h(1,:),'LineWidth',1,'LineStyle','-'); % Upper Whisker
        set(h(2,:),'LineWidth',1,'LineStyle','-'); % Lower Whisker
        set(h(3,:),'LineWidth',1,'LineStyle','-'); % Upper Adjacent Value
        set(h(4,:),'LineWidth',1,'LineStyle','-'); % Lower Adjacent Value
        set(h(5,:),'LineWidth',1,'LineStyle','-'); % Box
        set(h(6,:),'LineWidth',2,'LineStyle','-'); % Median


        hold on;

        % Overlay individual data points with colors matching cluster groups
        for cl = unique(grouping_data)'
            idx = (grouping_data == cl); % Find indices for each cluster
            sw{cl} = swarmchart(grouping_data(idx), plot_data(idx), 3, ...
                'filled', 'MarkerFaceColor', cluster_colors(cl, :), ...
                'MarkerFaceAlpha', 0.7, 'XJitter', 'density', 'XJitterWidth', 0.4);
        end

        T = table(chasingbycompetition_topleft_all(:), chasingbycompetition_topright_all(:), chasingbycompetition_bottomleft_all(:), chasingbycompetition_bottomright_all(:), ...
            'VariableNames', {'topleft', 'topright', 'bottomleft', 'bottomright'});
        writetable(T, fullfile(results_dir,'SourceData_ChasingByCompetition_AllDataPoints.csv'));

        %% 6) Plot quadrant proportions for Competition
        subplot(3,2,3)

        % Combine the four quadrant vectors into columns for notBoxPlot
        competitionMatrix = [competition_topleft, competition_topright, competition_bottomleft, competition_bottomright];
        nb=notBoxPlot(competitionMatrix);
        for nb_idx = 1:length(nb)
            nb(nb_idx).data.MarkerEdgeColor = 'none';
            nb(nb_idx).data.MarkerFaceColor = cluster_colors(nb_idx,:);
            nb(nb_idx).data.MarkerSize = 3;
            nb(nb_idx).mu.Color = cluster_colors(nb_idx,:);
            nb(nb_idx).sdPtch(1).FaceColor = cluster_colors(nb_idx,:);
            nb(nb_idx).sdPtch(1).FaceAlpha = 0.2;
            nb(nb_idx).semPtch(1).FaceColor = cluster_colors(nb_idx,:);
            nb(nb_idx).semPtch(1).FaceAlpha = 0.4;
            nb(nb_idx).sdPtch(1).EdgeColor = 'none';
            % nb(nb_idx).sdPtch(2).EdgeColor = 'none';
            nb(nb_idx).semPtch(1).EdgeColor = 'none';
            % nb(nb_idx).semPtch(2).EdgeColor = 'none';
        end

        % Improve axes
        ax2 = gca;
        ax2.XTick = 1:4;
        ax2.XTickLabel = {'Top-Left','Top-Right','Bottom-Left','Bottom-Right'};
        ylabel('Proportion of total competition');
        title('Competition (sorted by Competition DS)');
        ylim([0 1]);

        hold on; % To add permutation-test results on the figure



        %% 4) Plot quadrant proportions for Chasing (high-resolution data)
        subplot(3,2,4)

        % Combine the four quadrant vectors into columns for notBoxPlot
        competitionMatrix_all = [competition_topleft_all(:), competition_topright_all(:), competition_bottomleft_all(:), competition_bottomright_all(:)];
        % chasingbycompetitionMatrix_all = [chasingbycompetition_topleft_all(:).^(1/3), chasingbycompetition_topright_all(:).^(1/3), chasingbycompetition_bottomleft_all(:).^(1/3), chasingbycompetition_bottomright_all(:).^(1/3)];
        % % % t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        % % nexttile;
        cluster_colors = [0 0 255/2; 0 0 255/2; 0 0 255/2; 0 0 255/2; 0 0 255/2] / 255;
        plot_data = [];
        grouping_data = [];
        for q_idx = 1:size(competitionMatrix_all,2)
            plot_data = [plot_data;competitionMatrix_all(:,q_idx)];
            grouping_data = [grouping_data;q_idx*ones(size(competitionMatrix_all,1),1)];
        end
        h = boxplot(plot_data, grouping_data, 'Symbol', 'o', 'OutlierSize', 0.01, 'Whisker', 1, 'colors', cluster_colors, 'boxstyle', 'outline');
        % Find and modify the properties of different boxplot elements
        set(h(1,:),'LineWidth',1,'LineStyle','-'); % Upper Whisker
        set(h(2,:),'LineWidth',1,'LineStyle','-'); % Lower Whisker
        set(h(3,:),'LineWidth',1,'LineStyle','-'); % Upper Adjacent Value
        set(h(4,:),'LineWidth',1,'LineStyle','-'); % Lower Adjacent Value
        set(h(5,:),'LineWidth',1,'LineStyle','-'); % Box
        set(h(6,:),'LineWidth',2,'LineStyle','-'); % Median


        hold on;

        % Overlay individual data points with colors matching cluster groups
        for cl = unique(grouping_data)'
            idx = (grouping_data == cl); % Find indices for each cluster
            sw{cl} = swarmchart(grouping_data(idx), plot_data(idx), 3, ...
                'filled', 'MarkerFaceColor', cluster_colors(cl, :), ...
                'MarkerFaceAlpha', 0.7, 'XJitter', 'density', 'XJitterWidth', 0.4);
        end

        T = table(competition_topleft_all(:), competition_topright_all(:), competition_bottomleft_all(:), competition_bottomright_all(:), ...
            'VariableNames', {'topleft', 'topright', 'bottomleft', 'bottomright'});
        writetable(T, fullfile(results_dir,'SourceData_CompetitionByCompetition_AllDataPoints.csv'));

        %% 4) Plot quadrant proportions for Chasing (high-resolution data)
        subplot(3,2,5)

        % Combine the four quadrant vectors into columns for notBoxPlot
        chasingbychasingMatrix = [chasingbychasing_topleft, chasingbychasing_topright, chasingbychasing_bottomleft, chasingbychasing_bottomright];

        nb=notBoxPlot(chasingbychasingMatrix);
        for nb_idx = 1:length(nb)
            nb(nb_idx).data.MarkerEdgeColor = 'none';
            nb(nb_idx).data.MarkerFaceColor = cluster_colors(nb_idx,:);
            nb(nb_idx).data.MarkerSize = 3;
            nb(nb_idx).mu.Color = cluster_colors(nb_idx,:);
            nb(nb_idx).sdPtch(1).FaceColor = cluster_colors(nb_idx,:);
            nb(nb_idx).sdPtch(1).FaceAlpha = 0.2;
            nb(nb_idx).semPtch(1).FaceColor = cluster_colors(nb_idx,:);
            nb(nb_idx).semPtch(1).FaceAlpha = 0.4;
            nb(nb_idx).sdPtch(1).EdgeColor = 'none';
            % nb(nb_idx).sdPtch(2).EdgeColor = 'none';
            nb(nb_idx).semPtch(1).EdgeColor = 'none';
            % nb(nb_idx).semPtch(2).EdgeColor = 'none';
        end

        % Improve axes
        ax = gca;
        ax.XTick = 1:4;
        xlabels = {'Top-Left','Top-Right','Bottom-Left','Bottom-Right'};
        ax.XTickLabel = xlabels;
        ylabel('Proportion of total chasing');
        title('Chasing (sorted by Fraction of active chases)');
        ylim([0 1]);

        hold on; % To add permutation-test results on the figure

        %% 4) Plot quadrant proportions for Chasing (high-resolution data)
        subplot(3,2,6)

        % Combine the four quadrant vectors into columns for notBoxPlot
        chasingbychasingMatrix_all = [chasingbychasing_topleft_all(:), chasingbychasing_topright_all(:), chasingbychasing_bottomleft_all(:), chasingbychasing_bottomright_all(:)];
        % chasingbychasingMatrix_all = [chasingbychasing_topleft_all(:).^(1/3), chasingbychasing_topright_all(:).^(1/3), chasingbychasing_bottomleft_all(:).^(1/3), chasingbychasing_bottomright_all(:).^(1/3)];
        % % % t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        % % nexttile;
        cluster_colors = [0 0 255/2; 0 0 255/2; 0 0 255/2; 0 0 255/2; 0 0 255/2] / 255;
        plot_data = [];
        grouping_data = [];
        for q_idx = 1:size(chasingbychasingMatrix_all,2)
            plot_data = [plot_data;chasingbychasingMatrix_all(:,q_idx)];
            grouping_data = [grouping_data;q_idx*ones(size(chasingbychasingMatrix_all,1),1)];
        end
        h = boxplot(plot_data, grouping_data, 'Symbol', 'o', 'OutlierSize', 0.01, 'Whisker', 1, 'colors', cluster_colors, 'boxstyle', 'outline');
        % Find and modify the properties of different boxplot elements
        set(h(1,:),'LineWidth',1,'LineStyle','-'); % Upper Whisker
        set(h(2,:),'LineWidth',1,'LineStyle','-'); % Lower Whisker
        set(h(3,:),'LineWidth',1,'LineStyle','-'); % Upper Adjacent Value
        set(h(4,:),'LineWidth',1,'LineStyle','-'); % Lower Adjacent Value
        set(h(5,:),'LineWidth',1,'LineStyle','-'); % Box
        set(h(6,:),'LineWidth',2,'LineStyle','-'); % Median


        hold on;

        % Overlay individual data points with colors matching cluster groups
        for cl = unique(grouping_data)'
            idx = (grouping_data == cl); % Find indices for each cluster
            sw{cl} = swarmchart(grouping_data(idx), plot_data(idx), 3, ...
                'filled', 'MarkerFaceColor', cluster_colors(cl, :), ...
                'MarkerFaceAlpha', 0.7, 'XJitter', 'density', 'XJitterWidth', 0.4);
        end

        T = table(chasingbychasing_topleft_all(:), chasingbychasing_topright_all(:), chasingbychasing_bottomleft_all(:), chasingbychasing_bottomright_all(:), ...
            'VariableNames', {'topleft', 'topright', 'bottomleft', 'bottomright'});
        writetable(T, fullfile(results_dir,'SourceData_ChasingByChasing_AllDataPoints.csv'));


        % Save figures
        fig_filename_base = fullfile(results_dir, ['Quadrant_Analysis_SortedBy', chasing_sortingType, '_', day_ranges{day_idx}, '_', normalization]);
        [annot, srcInfo] = docDataSrc(fig,fullfile(results_dir),mfilename('fullpath'),logical(1));
        exportgraphics(fig, [fig_filename_base, '.pdf'], 'Resolution', 300);
        quadrant_quantification = array2table([chasingbycompetition_topleft,chasingbycompetition_topright,chasingbycompetition_bottomleft,chasingbycompetition_bottomright],'VariableNames',{'chasing_by_comp_TopLeft','chasing_by_comp_TopRight','chasing_by_comp_BottomLeft','chasing_by_comp_BottomRight'})
        quadrant_quantification = [cohortsTbl(:,1),quadrant_quantification];
        save(fullfile(results_dir, ['Quadrant_Quantification_ChasingSortedByComp_', day_ranges{day_idx}, '.mat']),'quadrant_quantification');
        chasing_direction = table(down_vs_up_chasing,down_vs_up_chasing_weighted,top_down_chasing,top_down_chasing_weighted,bottom_up_chasing,bottom_up_chasing_weighted,'VariableNames',{'down_vs_up_chasing','down_vs_up_chasing_weighted','top_down_chasing','top_down_chasing_weighted','bottom_up_chasing','bottom_up_chasing_weighted'});
        save(fullfile(results_dir, ['ChasingDirection_', day_ranges{day_idx}, '.mat']),'chasing_direction');
        quadrant_quantification_2 = array2table([competition_topleft,competition_topright,competition_bottomleft,competition_bottomright],'VariableNames',{'competition_by_comp_TopLeft','competition_by_comp_TopRight','competition_by_comp_BottomLeft','competition_by_comp_BottomRight'})
        quadrant_quantification_all = [quadrant_quantification, quadrant_quantification_2];
        writetable(quadrant_quantification_all,fullfile(results_dir,'SourceData_QuadrantQuantification.csv'));

        %% Plotting
        fig = figure('Visible', 'on');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.9]);

        myFields = fieldnames(myMatchMatrix);
        for subpl_idx = 1:length(myFields)
            subplot(2, 2, subpl_idx);
            imagesc(squeeze(nanmean(cat(3, myMatchMatrix.(myFields{subpl_idx}){:}), 3)));

            ax = gca;
            axis square;
            ax.XTick = 1:size(myMatchMatrix.(myFields{subpl_idx}){1}, 1);
            ax.YTick = 1:size(myMatchMatrix.(myFields{subpl_idx}){1}, 1);

            if contains(myFields{subpl_idx}, 'sortedByComp') && contains(myFields{subpl_idx}, 'Chasing')
                ax.XLabel.String = {'being chased', '(animals sorted by competition rank)'};
                ax.YLabel.String = {'chasings', '(animals sorted by competition rank)'};
            elseif contains(myFields{subpl_idx}, 'sortedByComp') && contains(myFields{subpl_idx}, 'Competition')
                ax.XLabel.String = {'losses', '(animals sorted by competition rank)'};
                ax.YLabel.String = {'wins', '(animals sorted by competition rank)'};
            elseif contains(myFields{subpl_idx}, 'sortedByChas') && contains(myFields{subpl_idx}, 'Chasing')
                ax.XLabel.String = {'being chased', '(animals sorted by chasing rank)'};
                ax.YLabel.String = {'chasings', '(animals sorted by chasing rank)'};
            elseif contains(myFields{subpl_idx}, 'sortedByChas') && contains(myFields{subpl_idx}, 'Competition')
                ax.XLabel.String = {'losses', '(animals sorted by chasing rank)'};
                ax.YLabel.String = {'wins', '(animals sorted by chasing rank)'};
            end

            % Title and colorbar
            title(strrep(myFields{subpl_idx}, '_', ' '), 'Interpreter', 'none');
            colorbar;
            ax.Colormap = parula;

            % Save source data
            writematrix(squeeze(nanmean(cat(3, myMatchMatrix.(myFields{subpl_idx}){:}), 3)), ...
                fullfile(results_dir, [myFields{subpl_idx}, '_', chasing_sortingType, '_', day_ranges{day_idx}, '_', normalization, '.csv']));
        end

        % Add a supertitle
        suptitle_str = sprintf('Sorting: %s, %s, Normalization: %s', chasing_sortingType, day_ranges{day_idx}, normalization);
        sgtitle(suptitle_str, 'Interpreter', 'none');

        % Save figures
        [annot, srcInfo] = docDataSrc(fig,fullfile(results_dir),mfilename('fullpath'),logical(1));
        fig_filename_base = fullfile(results_dir, ['Match_Matrices_SortedBy', chasing_sortingType, '_', day_ranges{day_idx}, '_', normalization]);
        exportgraphics(fig, [fig_filename_base, '.pdf'], 'Resolution', 300);
        exportgraphics(fig, [fig_filename_base, '.png'], 'Resolution', 300);
    end
end


%%
function myMatchMatrix = apply_normalization(myMatchMatrix, match_matrix_Competition, ...
    match_matrix_Chasing, sorting_idx_Comp, sorting_idx_Chas, cohort_idx, normalization)

% Empty NAN matrix
myMatchMatrix.Chasing_sortedByChas{cohort_idx}=nan(10,10);
myMatchMatrix.Chasing_sortedByComp{cohort_idx}=nan(10,10);
myMatchMatrix.Competition_sortedByChas{cohort_idx}=nan(10,10);
myMatchMatrix.Competition_sortedByComp{cohort_idx}=nan(10,10);

% Select normalization method
switch normalization
    case 'InteractionInDyad'
        normalize_func = @(mat, idx) mat(idx, idx) ./ (mat(idx, idx) + mat(idx, idx)');
    case 'maximum'
        normalize_func = @(mat, idx) mat(idx, idx) / max(max(mat(idx, idx)));
    case 'events'
        normalize_func = @(mat, idx) mat(idx, idx) / sum(sum(mat(idx, idx)));
    case 'none'
        normalize_func = @(mat, idx) mat(idx, idx);
end

% Apply normalization to sorted matrices
myMatchMatrix.Competition_sortedByComp{cohort_idx}(1:size(match_matrix_Competition,1),1:size(match_matrix_Competition,2)) = normalize_func(match_matrix_Competition, sorting_idx_Comp);
myMatchMatrix.Competition_sortedByChas{cohort_idx}(1:size(match_matrix_Competition,1),1:size(match_matrix_Competition,2)) = normalize_func(match_matrix_Competition, sorting_idx_Chas);
myMatchMatrix.Chasing_sortedByComp{cohort_idx}(1:size(match_matrix_Chasing,1),1:size(match_matrix_Chasing,2)) = normalize_func(match_matrix_Chasing, sorting_idx_Comp);
myMatchMatrix.Chasing_sortedByChas{cohort_idx}(1:size(match_matrix_Chasing,1),1:size(match_matrix_Chasing,2)) = normalize_func(match_matrix_Chasing, sorting_idx_Chas);
end
