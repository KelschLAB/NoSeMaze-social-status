function [Results, lme_model] = permutation_and_lme_analysis(yData, group, animal_ID, varargin)
% Perform permutation by swapping entire animal trial sets and GLM on group differences.
% Supports 'mean', 'median', and 'ranks' modes.
%
% Inputs:
%   yData - Response variable
%   group - Categorical group (WT vs KO)
%   animal_ID - Animal identifiers
% Optional:
%   'n_perm' - Number of permutations (default 10000)
%   'method' - 'mean', 'median' or 'ranks' (default 'mean')
%   'run_lme' - Logical, run GLM as control (default true)
%
% Output:
%   Results - Structure with p-values, effect sizes, and GLM if requested

%% Parse Inputs
p = inputParser;
addOptional(p, 'n_perm', 10000);
addOptional(p, 'method', 'mean');
addOptional(p, 'run_lme', true);
parse(p, varargin{:});
n_perm = p.Results.n_perm;
method = p.Results.method;
run_lme = p.Results.run_lme;

%% Prepare Data
yData = yData(:);
group = group(:);
animal_ID = animal_ID(:);

% Ensure string format
if iscell(animal_ID)
    animal_ID = string(animal_ID);
end
if iscell(group)
    group = string(group);
end

% Prepare numeric group labels
group_numeric = double(group == 'OxtKO'); % 0 = WT, 1 = KO

% Group indices
idx_WT = group_numeric == 0;
idx_KO = group_numeric == 1;

%% Effect Size Calculations
% Mean and Median differences
mean_diff = mean(yData(idx_KO)) - mean(yData(idx_WT));
median_diff = median(yData(idx_KO)) - median(yData(idx_WT));

% Cohen's d
n1 = sum(idx_WT);
n2 = sum(idx_KO);
s1 = std(yData(idx_WT));
s2 = std(yData(idx_KO));
pooled_sd = sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2));
cohen_d = mean_diff / pooled_sd;

% Cliff's delta
delta = cliffs_delta(yData(idx_WT), yData(idx_KO));

%% Observed difference based on method
disp('Calculating observed difference...');

switch method
    case 'mean'
        diff_obs = mean_diff;
        effect_size = cohen_d;
    case 'median'
        diff_obs = median_diff;
        effect_size = delta;  % optional: can add Hodges–Lehmann if desired
    case 'ranks'
        diff_obs = cliffs_delta(yData(idx_WT), yData(idx_KO));
        effect_size = diff_obs;  % delta is the effect size
    otherwise
        error('Invalid method. Choose ''mean'', ''median'', or ''ranks''.');
end

%% Skewness Calculations
% Skewness of the original yData for each group (WT vs KO)
skewness_WT = skewness(yData(idx_WT));
skewness_KO = skewness(yData(idx_KO));
disp(['Skewness of WT group: ', num2str(skewness_WT)]);
disp(['Skewness of KO group: ', num2str(skewness_KO)]);

% Skewness for log-transformed data (if applicable)
if run_lme
    shifted_y = yData - min(yData) + 1;  % Log shift for positive values
    log_yData = log(shifted_y);
    skewness_log = skewness(log_yData);
    disp(['Skewness of log-transformed yData: ', num2str(skewness_log)]);
end

%% Build animal mapping
disp('Preparing permutation structure...');
animals = unique(animal_ID);
animal_groups = zeros(length(animals),1);

for i = 1:length(animals)
    idx = find(animal_ID == animals(i));
    animal_groups(i) = group_numeric(idx(1)); % All trials have same group
end

%% Permutation: shuffle animal groups
disp('Running block-level permutation...');
perm_diffs = zeros(n_perm,1);

for perm_idx = 1:n_perm
    shuffled_labels = animal_groups(randperm(length(animal_groups))); % Shuffle animal labels

    % Assign new group labels back to trials
    permuted_group_numeric = zeros(size(group_numeric));
    for i = 1:length(animals)
        permuted_group_numeric(animal_ID == animals(i)) = shuffled_labels(i);
    end

    % Recalculate difference
    switch method
        case 'mean'
            perm_diffs(perm_idx) = nanmean(yData(permuted_group_numeric==1)) - ...
                                   nanmean(yData(permuted_group_numeric==0));
        case 'median'
            perm_diffs(perm_idx) = nanmedian(yData(permuted_group_numeric==1)) - ...
                                   nanmedian(yData(permuted_group_numeric==0));
        case 'ranks'
            perm_diffs(perm_idx) = cliffs_delta(yData(permuted_group_numeric==0), ...
                                        yData(permuted_group_numeric==1));
    end
end

% Skewness for permutation differences
skewness_perm = skewness(perm_diffs);
disp(['Skewness of permutation distribution: ', num2str(skewness_perm)]);

% Compute permutation p-value (two-tailed)
pval_permutation = mean(abs(perm_diffs) >= abs(diff_obs));

%% GLM on full trial data
if run_lme
    disp('Running GLM...');
    % Raw model
    lme_tbl = table(group, animal_ID, yData, ...
        'VariableNames', {'Group','animal_ID','Response'});
    lme_tbl.Group = categorical(lme_tbl.Group);
    lme_tbl.animal_ID = categorical(lme_tbl.animal_ID);
    lme_model = fitlme(lme_tbl, 'Response ~ Group + (Group|animal_ID)');
    lme_coeff = lme_model.Coefficients;
    pval_lme = lme_coeff.pValue(2); % Group effect p-value

    % Log-transformed model
    % shifted_y = yData - min(yData) + 1;
    % log_yData = log(shifted_y);
    shifted_y = yData - min(yData) + 0.001;
    log_yData = boxcox(shifted_y);%figure; histogram(log_yData);
    lme_log_tbl = table(group, animal_ID, log_yData, ...
        'VariableNames', {'Group','animal_ID','Response'});
    lme_log_tbl.Group = categorical(lme_log_tbl.Group);
    lme_log_tbl.animal_ID = categorical(lme_log_tbl.animal_ID);
    log_lme_model = fitlme(lme_log_tbl, 'Response ~ Group + (Group|animal_ID)');
    log_lme_coeff = log_lme_model.Coefficients;
    pval_lme_log = log_lme_coeff.pValue(2);
else
    lme_model = [];
    log_lme_model = [];
    pval_lme = NaN;
    pval_lme_log = NaN;
end


%% Package Results
Results.diff_obs = diff_obs;
Results.pval_permutation = pval_permutation;
Results.mean_diff = mean_diff;
Results.median_diff = median_diff;
Results.cohen_d = cohen_d;
Results.cliffs_delta = delta;
Results.effect_size = effect_size;
Results.method_used = method;
Results.skewness_WT = skewness_WT;
Results.skewness_KO = skewness_KO;
Results.skewness_perm = skewness_perm;
if run_lme
    Results.pval_lme = pval_lme;
    Results.pval_lme_log = pval_lme_log;
    Results.lme_model = lme_model;
    Results.log_lme_model = log_lme_model;
    Results.skewness_log = skewness_log;
end


end

%% Cliff's delta helper function
function delta = cliffs_delta(x, y)
% Compute Cliff's delta: difference in probability dominance
n_x = numel(x);
n_y = numel(y);
[xx, yy] = ndgrid(x, y);
delta = (sum(xx(:) < yy(:)) - sum(xx(:) > yy(:))) / (n_x * n_y);
end
