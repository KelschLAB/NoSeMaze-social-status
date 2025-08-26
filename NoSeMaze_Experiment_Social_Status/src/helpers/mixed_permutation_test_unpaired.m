function [p_value, diff_obs] = mixed_permutation_test_unpaired(xData, yData, animal_id_x, animal_id_y, n_perm, global_swap_ratio, method, mode, run_glm_control)
% xData: Data for cluster 1
% yData: Data for cluster 2
% animal_id_x: Animal IDs corresponding to xData
% animal_id_y: Animal IDs corresponding to yData
% n_perm: Number of permutations
% global_swap_ratio: Proportion of trials swapped (only used if mode == 'trial')
% method: 'mean' or 'median'
% mode: 'trial' or 'animal'

if iscell(animal_id_x), animal_id_x = string(animal_id_x); end
if iscell(animal_id_y), animal_id_y = string(animal_id_y); end

all_data = [xData; yData];
all_animal_id = [animal_id_x; animal_id_y];
all_group_labels = [ones(size(xData)); 2 * ones(size(yData))];

if strcmp(mode, 'animal')
    unique_animals = unique(all_animal_id);

    agg_data = [];
    agg_labels = [];

    rng('shuffle') % Optional: random seed for random assignments

    for a = unique_animals'
        a = string(a);
        x_vals = xData(animal_id_x == a);
        y_vals = yData(animal_id_y == a);

        if isempty(x_vals) && ~isempty(y_vals)
            % Only y-data
            if strcmp(method, 'mean')
                agg_data(end+1,1) = mean(y_vals);
            else
                agg_data(end+1,1) = median(y_vals);
            end
            agg_labels(end+1,1) = 2;
        elseif isempty(y_vals) && ~isempty(x_vals)
            % Only x-data
            if strcmp(method, 'mean')
                agg_data(end+1,1) = mean(x_vals);
            else
                agg_data(end+1,1) = median(x_vals);
            end
            agg_labels(end+1,1) = 1;
        elseif ~isempty(x_vals) && ~isempty(y_vals)
            % Animal has both -> randomly assign to one group
            if rand() > 0.5
                % Assign to x-group
                if strcmp(method, 'mean')
                    agg_data(end+1,1) = mean(x_vals);
                else
                    agg_data(end+1,1) = median(x_vals);
                end
                agg_labels(end+1,1) = 1;
            else
                % Assign to y-group
                if strcmp(method, 'mean')
                    agg_data(end+1,1) = mean(y_vals);
                else
                    agg_data(end+1,1) = median(y_vals);
                end
                agg_labels(end+1,1) = 2;
            end
        end
    end

    % Observed difference
    diff_obs = mean(agg_data(agg_labels==1)) - mean(agg_data(agg_labels==2));

    % Permutation
    perm_diffs = zeros(n_perm,1);
    for i = 1:n_perm
        perm_labels = agg_labels(randperm(length(agg_labels)));
        perm_diffs(i) = mean(agg_data(perm_labels==1)) - mean(agg_data(perm_labels==2));
    end

elseif strcmp(mode, 'trial')
    % Trial-level permutation with clustered swap
    all_animals = union(unique(animal_id_x), unique(animal_id_y));

    % Compute observed difference
    if strcmp(method, 'mean')
        diff_obs = mean(xData) - mean(yData);
    elseif strcmp(method, 'median')
        diff_obs = median(xData) - median(yData);
    else
        error('Invalid method. Choose either ''mean'' or ''median''.');
    end

    perm_diffs = zeros(n_perm, 1);
    for i = 1:n_perm
        perm_labels = all_group_labels; % Start with original labels

        % Shuffle within animals
        for a = all_animals'
            idx = find(all_animal_id == a);
            if ~isempty(idx) && length(unique(perm_labels(idx))) > 1
                swap_idx = rand(size(idx)) > 0.5;
                if any(swap_idx)
                    temp = perm_labels(idx(swap_idx));
                    perm_labels(idx(swap_idx)) = 3 - temp;
                end
            end
        end

        % Global swap
        num_global_swaps = round(length(all_data) * global_swap_ratio);
        global_swap_idx = randperm(length(all_data), num_global_swaps);
        temp = perm_labels(global_swap_idx);
        perm_labels(global_swap_idx) = 3 - temp;

        % Compute permuted difference
        perm_x = all_data(perm_labels == 1);
        perm_y = all_data(perm_labels == 2);

        if strcmp(method, 'mean')
            perm_diffs(i) = mean(perm_x) - mean(perm_y);
        elseif strcmp(method, 'median')
            perm_diffs(i) = median(perm_x) - median(perm_y);
        end
    end
else
    error('Invalid mode. Choose either ''trial'' or ''animal''.');
end
% Two-tailed p-value
p_value = mean(abs(perm_diffs) >= abs(diff_obs));

% --- Now run optional GLM control ---
if exist('run_glm_control', 'var') && run_glm_control
    % Prepare data for GLM
    glm_data = all_data;
    glm_group = categorical(all_group_labels); % Group 1 vs 2

    % Optional: transform if necessary
    % glm_data = log(glm_data ./ (1 - glm_data)); % Only if you want logit

    % Fit GLM
    glm_model = fitglm(glm_group, glm_data, 'Distribution', 'normal', 'Link', 'identity'); % Simple linear model
    glm_table = glm_model.Coefficients;
    pval_glm = glm_table.pValue(2); % p-value for Group effect

    disp(['GLM p-value for group effect: ', num2str(pval_glm)]);
    disp(['Permutation p-value for group effect: ', num2str(p_value)]);
end

% End function properly
end


%
%
% function [p_value,diff_obs] = mixed_permutation_test_unpaired(xData, yData, animal_id_x, animal_id_y, n_perm, global_swap_ratio, method)
%     % xData: Data for cluster 1 (e.g., CS+ latency in cluster 1)
%     % yData: Data for cluster 2 (e.g., CS+ latency in cluster 2)
%     % animal_id_x: Animal IDs corresponding to xData
%     % animal_id_y: Animal IDs corresponding to yData
%     % n_perm: Number of permutations
%     % global_swap_ratio: Proportion of trials swapped across animals (0 to 1)
%     % method: 'mean' or 'median' to determine which statistic to use
%
%     % Convert animal IDs from cell arrays to strings or categorical arrays
%     if iscell(animal_id_x), animal_id_x = string(animal_id_x); end
%     if iscell(animal_id_y), animal_id_y = string(animal_id_y); end
%
%     % Step 1: Find all unique animals across both groups (not just matched ones)
%     all_animals = union(unique(animal_id_x), unique(animal_id_y));
%
%     % Step 2: Combine all data from both clusters
%     all_data = [xData; yData];
%     all_animal_id = [animal_id_x; animal_id_y];
%     all_group_labels = [ones(size(xData)); 2 * ones(size(yData))]; % 1 = Cluster 1, 2 = Cluster 2
%
%     % Step 3: Compute observed difference using specified method
%     if strcmp(method, 'mean')
%         diff_obs = mean(xData) - mean(yData);
%     elseif strcmp(method, 'median')
%         diff_obs = median(xData) - median(yData);
%     else
%         error('Invalid method. Choose either ''mean'' or ''median''.');
%     end
%
%     perm_diffs = zeros(n_perm, 1);
%
%     for i = 1:n_perm
%         perm_labels = all_group_labels;  % Copy original labels
%
%         % Step 4: Shuffle within animals (for animals appearing in both groups)
%         for a = all_animals'
%             idx = find(all_animal_id == a);
%             if ~isempty(idx) && length(unique(perm_labels(idx))) > 1  % Ensure at least 2 conditions exist
%                 swap_idx = rand(size(idx)) > 0.5; % Randomly swap 50% of trials
%                 if any(swap_idx)  % Only swap if there are valid trials to swap
%                     temp = perm_labels(idx(swap_idx));
%                     perm_labels(idx(swap_idx)) = 3 - temp; % Flip 1↔2
%                 end
%             end
%         end
%
%         % Step 5: Global shuffle (swapping between all animals, even unpaired)
%         num_global_swaps = round(length(all_data) * global_swap_ratio); % Variable swap percentage
%         global_swap_idx = randperm(length(all_data), num_global_swaps);
%         temp = perm_labels(global_swap_idx);
%         perm_labels(global_swap_idx) = 3 - temp; % Flip 1↔2
%
%         % Step 6: Compute permuted statistic using chosen method
%         perm_x = all_data(perm_labels == 1);
%         perm_y = all_data(perm_labels == 2);
%
%         if strcmp(method, 'mean')
%             perm_diffs(i) = mean(perm_x) - mean(perm_y);
%         elseif strcmp(method, 'median')
%             perm_diffs(i) = median(perm_x) - median(perm_y);
%         end
%     end
%
%     % Step 7: Compute p-value (two-tailed test)
%     p_value = mean(abs(perm_diffs) >= abs(diff_obs));
% end
