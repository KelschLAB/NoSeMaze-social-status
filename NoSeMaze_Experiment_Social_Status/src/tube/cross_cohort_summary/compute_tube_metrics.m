function metrics = compute_tube_metrics(match_matrix, plot_dir)
% Computes key metrics (DSz, rank, wins/losses) from a match matrix.
%
% Parameters:
% match_matrix - Square matrix of wins/losses between animals
%
% Returns:
% metrics - Struct containing calculated metrics

metrics = struct();

% Number of wins (rows) and losses (columns)
metrics.N_Wins = sum(match_matrix, 2); % Sum of rows (wins)
metrics.N_Losses = sum(match_matrix, 1)'; % Sum of columns (losses)
metrics.N_Events = metrics.N_Wins + metrics.N_Losses; % Sum of wins and losses
metrics.N_WinsMinusLosses_normalized = (metrics.N_Wins - metrics.N_Losses)./(metrics.N_Wins + metrics.N_Losses); % normalized

% Fractions of wins and losses
total_matches = sum(match_matrix(:));
metrics.Fraction_Wins = metrics.N_Wins ./ total_matches; % Fraction of wins
metrics.Fraction_Losses = metrics.N_Losses ./ total_matches; % Fraction of losses

% Compute David's Score (DS) and its z-score (DSz)
DS_info = compute_DS_from_match_matrix(match_matrix); % Pre-existing helper function
metrics.DSz = zscore(DS_info.DS)'; % Standardized David's Score (z-scored)
metrics.DS = (DS_info.DS)'; % Standardized David's Score (z-scored)

% Compute rank based on DSz
[~, sorted_indices] = sort(metrics.DSz, 'descend'); % Sort DSz in descending order
metrics.Rank = zeros(size(metrics.DSz)); % Initialize rank array
metrics.Rank(sorted_indices) = 1:length(metrics.DSz); % Assign ranks based on sorting

%% Tranformation of metrics (searching for optimal lambda and transform
% data)
transformation_metrics = {'N_Wins','N_Losses','N_Events','Fraction_Wins','Fraction_Losses'};
% Loop over transformation_metrics
for ii = 1:length(transformation_metrics)
    % data definition
    data = metrics.(transformation_metrics{ii});
    % Add a small constant to avoid zero or negative values
    if min(data)<0
        epsilon = ceil(abs(min(data)));
    else
        epsilon = 1e-3; % Small constant
    end
    data_adj = data + epsilon;

    % Define a range of lambda values to test
    lambda_values = -10:0.01:10;
    % Initialize variable to store the maximum log-likelihood and the best lambda
    max_log_likelihood = -Inf;
    best_lambda = NaN;
    % Iterate over each lambda value
    for lambda = lambda_values
        if lambda == 0
            % For lambda = 0, use log transformation
            transformed_data = log(data_adj);
        else
            % Box-Cox transformation for lambda != 0
            transformed_data = (data_adj .^ lambda - 1) / lambda;
        end
        % Calculate the log-likelihood
        n = length(data_adj);
        log_likelihood = -n/2 * log(var(transformed_data)) + (lambda - 1) * sum(log(data_adj));
        % Update the best lambda if the current log-likelihood is higher
        if log_likelihood > max_log_likelihood
            max_log_likelihood = log_likelihood;
            best_lambda = lambda;
        end
    end
    % Display the optimal lambda
    disp('Optimal Lambda:');
    disp(best_lambda);

    if best_lambda == 0
        transformed_data_boxcox = data_adj;
        transformed_data_log10 = log10(data);
        transformed_data_cuberoot = data .^ (1/3);
    else
        transformed_data_boxcox = (data_adj .^ best_lambda - 1) / best_lambda;
        transformed_data_log10 = log10(data);
        transformed_data_cuberoot = data .^ (1/3);
        % if boxcox is not working, all values might be equal, in this
        % case, use the standard data
        if isscalar(unique(transformed_data_boxcox))
            transformed_data_boxcox = data_adj;
        end
    end
    % 2.2 test for normality with transformed data
    [data_transformation(ii).H_sw_transformed_boxcox, data_transformation(ii).pValue_sw_transformed_boxcox, ~] = swtest(transformed_data_boxcox);
    [data_transformation(ii).H_ks_transformed_boxcox, data_transformation(ii).pValue_ks_transformed_boxcox, ~, ~] = kstest(transformed_data_boxcox);

    [data_transformation(ii).H_sw_log10, data_transformation(ii).pValue_sw_log10, ~] = swtest(transformed_data_log10);
    [data_transformation(ii).H_ks_log10, data_transformation(ii).pValue_ks_log10, ~, ~] = kstest(transformed_data_log10);

    [data_transformation(ii).H_sw_cuberoot, data_transformation(ii).pValue_sw_cuberoot, ~] = swtest(transformed_data_cuberoot);
    [data_transformation(ii).H_ks_cuberoot, data_transformation(ii).pValue_ks_cuberoot, ~, ~] = kstest(transformed_data_cuberoot);

    [data_transformation(ii).skewness_original] = skewness(data_adj);
    [data_transformation(ii).skewness_cuberoot] = skewness(transformed_data_cuberoot);
    [data_transformation(ii).skewness_log10] = skewness(transformed_data_log10);
    [data_transformation(ii).skewness_boxcox] = skewness(transformed_data_boxcox);

    data_transformation(ii).transformed_data_boxcox = transformed_data_boxcox;
    data_transformation(ii).transformed_data_log10 = transformed_data_log10;
    data_transformation(ii).transformed_data_cuberoot = transformed_data_cuberoot;
    data_transformation(ii).best_lambda = best_lambda;

    % 2.3 redo plotting
    f=figure('Visible','off');
    histogram(data_transformation(ii).transformed_data_boxcox,'EdgeColor','none');
    box off
    title([transformation_metrics{ii} ' (transformed boxcox)'],'Interpreter','none');
    saveas(f,fullfile(plot_dir,[transformation_metrics{ii},'_transformedBoxCox.png']),'png');
    saveas(f,fullfile(plot_dir,[transformation_metrics{ii},'_transformedBoxCox.pdf']),'pdf');
    close all;

    f=figure('Visible','off');
    histogram(data_transformation(ii).transformed_data_log10,'EdgeColor','none');
    box off
    title([transformation_metrics{ii} ' (transformed log10)'],'Interpreter','none');
    saveas(f,fullfile(plot_dir,[transformation_metrics{ii},'_transformedLog10.png']),'png');
    saveas(f,fullfile(plot_dir,[transformation_metrics{ii},'_transformedLog10.pdf']),'pdf');
    close all;

    f=figure('Visible','off');
    histogram(data_transformation(ii).transformed_data_cuberoot,'EdgeColor','none');
    box off
    title([transformation_metrics{ii} ' (transformed cuberoot)'],'Interpreter','none');
    saveas(f,fullfile(plot_dir,[transformation_metrics{ii},'_transformedcuberoot.png']),'png');
    saveas(f,fullfile(plot_dir,[transformation_metrics{ii},'_transformedcuberoot.pdf']),'pdf');
    close all;

    % preparation
    metrics.(['log10_' transformation_metrics{ii}]) = data_transformation(ii).transformed_data_log10;
    metrics.(['boxcox_' transformation_metrics{ii}]) = data_transformation(ii).transformed_data_boxcox;
    metrics.(['cuberoot_' transformation_metrics{ii}]) = data_transformation(ii).transformed_data_cuberoot ;

end
end
