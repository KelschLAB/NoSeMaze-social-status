%% analyze_time_to_criterion.m
% Author: Jonathan Reinwald
% Date: March 2025

%% Pre-Clearing
clear; clc; close all;

%% Set Directories
main_dir = 'myRootPath/NoSeMaze_Experiment_Social_Status'; % --> NOTE: Replace with own root directory
processed_dir = fullfile(main_dir, 'data', 'processed');
results_dir = fullfile(main_dir, 'results', 'figures', 'cross_cohort', 'lick_performance');
if ~isfolder(results_dir)
    mkdir(results_dir);
end
summary_file = fullfile(processed_dir, 'cross_cohort_files', 'summary_data.mat');

% Add required paths
addpath(genpath(fullfile(main_dir, 'src', 'helpers')));
% addpath(genpath(fullfile(main_dir, 'src', 'analysis', 'associations_between_metrics')));

%% Load Summary Data
if ~isfile(summary_file)
    error('Summary file not found. Please run `main_build_cross_cohort_summary` first.');
end
load(summary_file, 'summary_data');

%% Cohort Selection
cohortsTbl = readtable(fullfile(main_dir, 'config', 'cohorts_info.csv'));

%% Define Parameters
day_ranges = {'D1_End'}; % Select day_ranges of interest
cuberoot_transformation = 0;

% select cohorts
cohortsTbl = cohortsTbl(cohortsTbl.use_lickport == 1, :);

% Filter input data (cohort selection)
fields = fieldnames(summary_data.lickport);
for field_idx = 1:numel(fields)
    summary_data.lickport.(fields{field_idx}) = summary_data.lickport.(fields{field_idx})(ismember(summary_data.lickport.(fields{field_idx}).cohort,cohortsTbl.cohort),:);
end

%% Iterate Through day_ranges
for day_range_idx = 1:length(day_ranges)

    % clearing
    clear day_range data

    % Current day range
    day_range = day_ranges{day_range_idx};

    % Current data
    data = summary_data.lickport.(day_range);
    % Only consider the first repetition
    data = data(data.repetition==1,:);

    % Data curation
    phase_selection = 1:6;
    data_matrix_csplus = [];
    data_matrix_csminus = [];
    for phase_idx = 1:length(phase_selection)
        if cuberoot_transformation == 1
            data_matrix_csplus = [data_matrix_csplus,data.(['cuberoot_cs_plus_time_to_criterion_phase' num2str(phase_selection(phase_idx))])];
            data_matrix_csminus = [data_matrix_csminus,data.(['cuberoot_cs_minus_time_to_criterion_phase' num2str(phase_selection(phase_idx))])];
        else
            data_matrix_csplus = [data_matrix_csplus,data.(['cs_plus_time_to_criterion_phase' num2str(phase_selection(phase_idx))])];
            data_matrix_csminus = [data_matrix_csminus,data.(['cs_minus_time_to_criterion_phase' num2str(phase_selection(phase_idx))])];
        end
        phase_names{phase_idx} = ['phase_' num2str(phase_selection(phase_idx))];
    end
    % figure definition
    fig = figure('visible','on');
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.6]);
    % loop over subplots
    for subplot_idx = 1:2
        hold on;
        % plot data
        if subplot_idx==1
            plot_data = data_matrix_csplus;
        elseif subplot_idx==2
            plot_data = data_matrix_csminus;
        end
        % shadedErrorBar
        sd = shadedErrorBar(1:size(plot_data,2),nanmean(plot_data),SEM_calc(plot_data));
        sd.edge(1).Color='none';
        sd.edge(2).Color='none';
        if subplot_idx==1
            sd.mainLine.Color=[0.1 0.1 0.1];
            sd.patch.FaceColor=[0.45 0.45 0.45];
        else
            sd.mainLine.Color=[1 0 0];
            sd.patch.FaceColor=[0.45 0 0];
        end
        % % % subplot
        % % % subplot(1,2,subplot_idx);
        % % % notBoxPlot
        % % % nb = notBoxPlot(plot_data,'jitter',0.5);
        % % % color definitions
        % % % for ib = 1:length(nb)
        % % %     nb(ib).data.MarkerSize=6;
        % % %     nb(ib).data.MarkerEdgeColor='none';
        % % %     nb(ib).semPtch.EdgeColor='none';
        % % %     nb(ib).sdPtch.EdgeColor='none';
        % % %     nb(ib).data.MarkerFaceColor= [0.3 0.3 0.3];
        % % %     nb(ib).mu.Color= [0.1 0.1 0.1];
        % % %     nb(ib).semPtch.FaceColor= [0.45 0.45 0.45];
        % % %     nb(ib).sdPtch.FaceColor= [0.6 0.6 0.6];
        % % % end
        % Customize axis
        box('off');
        ax=gca;
        if cuberoot_transformation == 1
            ax.YLim=[0,12];
        else
            ax.YLim=[0,300];
        end
        % Customize plot
        ax.XLabel.String = 'phase';
        ax.XLabel.Interpreter = 'none';
        if cuberoot_transformation == 1
            ax.YLabel.String = {'trials until 80% criterion','(cuberoot transformed)'};
        else
            ax.YLabel.String = 'trials until 80% criterion';
        end
        ax.YLabel.Interpreter = 'none';
        ax.FontSize = 20;
        % title
        if subplot_idx==2
            title('CS, repetition 1','Interpreter','none');
            legend({'CS+','CS-'})
        end

        % Source data
        SourceDataTable = array2table(plot_data,'VariableNames',phase_names);
        if cuberoot_transformation == 1
            if subplot_idx==1
                writetable(SourceDataTable,fullfile(results_dir, 'NotBoxPlot_TimeToCriterion_cuberoot_CSplus.csv'));
            elseif subplot_idx==2
                writetable(SourceDataTable,fullfile(results_dir, 'NotBoxPlot_TimeToCriterion_cuberoot_CSminus.csv'));
            end
        else
            if subplot_idx==1
                writetable(SourceDataTable,fullfile(results_dir, 'NotBoxPlot_TimeToCriterion_CSplus.csv'));
            elseif subplot_idx==2
                writetable(SourceDataTable,fullfile(results_dir, 'NotBoxPlot_TimeToCriterion_CSminus.csv'));
            end
        end
    end
    % Save Plot
    [annot, srcInfo] = docDataSrc(fig,fullfile(results_dir),mfilename('fullpath'),logical(1));
    if cuberoot_transformation == 1
        exportgraphics(fig, fullfile(results_dir, 'NotBoxPlot_TimeToCriterion_cuberoot.pdf'));
    else
        exportgraphics(fig, fullfile(results_dir, 'NotBoxPlot_TimeToCriterion.pdf'));
    end
    close(fig);
end

