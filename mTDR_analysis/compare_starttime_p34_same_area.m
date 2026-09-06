%% compare_starttime_p34_same_area.m
% Compare onset/start time between two variables p=3 and p=4 within the same brain area.
%
% This script is modified from the original MST-vs-LIP comparison script.
% Original logic:
%   same p, compare MST_starttime vs LIP_starttime
%
% New logic:
%   same brain area, compare p3_starttime vs p4_starttime
%
% Statistical tests:
% 1) Bootstrap difference distribution
%       diff_boot = p4_boot - p3_boot
% 2) Paired label-swap permutation test
%       labels p3 and p4 are randomly swapped within each paired bootstrap repeat
%
% Interpretation of diff_boot:
%   diff_boot > 0: p4 starts later than p3
%   diff_boot < 0: p4 starts earlier than p3
%
% Required input files, for example:
%   ./timedata/Fp3MSTtime.mat, containing MST_starttime
%   ./timedata/Fp4MSTtime.mat, containing MST_starttime
%   ./timedata/Fp3LIPtime.mat, containing LIP_starttime
%   ./timedata/Fp4LIPtime.mat, containing LIP_starttime

clear; clc;

%% -----------------------------
%  Parameters
%  -----------------------------
monkey = 'D';
pair_p = [3, 4];          % compare p=3 vs p=4
brain_area_list = {'MST','LIP'};  % use {'MST'} or {'LIP'} if only one area is needed
nperm = 5000;
alpha = 0.05;

% x-axis range for onset/start time plot
xlim_range = [200, 1300];

% whether to save results
save_result = 0;%true;

%% -----------------------------
%  Color settings
%  -----------------------------
color_{1,1} = [205,0,0;0,0,0;0,100,0]/255;          % heading
color_{2,1} = [0,191,255;155,48,255]/255;           % color
color_{3,1} = [139,10,80;72,61,139]/255;            % choice / perceptual choice
color_{4,1} = [0,0,128;205,133,0]/255;              % saccade / saccadic choice
color_{5,1} = [0,0,139;205,133,0]/255;              % last saccade
color_{6,1} = [0,0,1;0,0,0.9;0,0,0.8;0,0,0.7;0,0,0.6]; % coherence
color_{7,1} = [205,0,0;0,0,0;0,100,0]/255;          % heading

p3_color = color_{3,1}(1,:);
p4_color = [30,114,255]/255;

all_results = struct();

%% -----------------------------
%  Main loop across brain areas
%  -----------------------------
for iarea = 1:numel(brain_area_list)

    brain_area = brain_area_list{iarea};

    fprintf('\n\n========================================\n');
    fprintf('Compare p%d vs p%d within %s, monkey %s\n', pair_p(1), pair_p(2), brain_area, monkey);
    fprintf('========================================\n');

    %% -----------------------------
    %  Load saved start-time files
    %  -----------------------------
    file_p1 = ['./timedata/', monkey, 'p', num2str(pair_p(1)), brain_area, 'time.mat'];
    file_p2 = ['./timedata/', monkey, 'p', num2str(pair_p(2)), brain_area, 'time.mat'];

    starttime_p1 = local_load_starttime(file_p1, brain_area);
    starttime_p2 = local_load_starttime(file_p2, brain_area);

    %% -----------------------------
    %  Make sure data are cell arrays
    %  -----------------------------
    if ~iscell(starttime_p1)
        starttime_p1 = num2cell(starttime_p1);
    end

    if ~iscell(starttime_p2)
        starttime_p2 = num2cell(starttime_p2);
    end

    %% -----------------------------
    %  Remove empty bootstrap samples
    %  -----------------------------
    valid_p1 = ~cellfun(@isempty, starttime_p1);
    valid_p2 = ~cellfun(@isempty, starttime_p2);

    % Paired bootstrap comparison:
    % only keep repeats where both p3 and p4 have valid onset times.
    valid_pair = valid_p1 & valid_p2;

    p1_boot = cell2mat(starttime_p1(valid_pair));
    p2_boot = cell2mat(starttime_p2(valid_pair));

    p1_boot = p1_boot(:);
    p2_boot = p2_boot(:);

    B = length(p1_boot);

    if B < 5
        error('Too few valid paired bootstrap samples for %s: B = %d', brain_area, B);
    end

    if length(p1_boot) ~= length(p2_boot)
        error('p1_boot and p2_boot have different lengths.');
    end

    fprintf('\nNumber of valid paired bootstrap samples = %d\n', B);

    %% -----------------------------
    %  Descriptive statistics for each variable p
    %  -----------------------------
    p1_mean = mean(p1_boot);
    p2_mean = mean(p2_boot);

    p1_median = median(p1_boot);
    p2_median = median(p2_boot);

    % Since p1_boot and p2_boot are already bootstrap distributions,
    % these percentile ranges describe bootstrap uncertainty of each variable's onset.
    p1_boot_CI = prctile(p1_boot, [100*alpha/2, 100*(1-alpha/2)]);
    p2_boot_CI = prctile(p2_boot, [100*alpha/2, 100*(1-alpha/2)]);

    %% -----------------------------
    %  Bootstrap difference distribution
    %  -----------------------------
    % Positive value:
    %   p4 starts later than p3
    % Negative value:
    %   p4 starts earlier than p3
    diff_boot = p2_boot - p1_boot;

    diff_mean = mean(diff_boot);
    diff_median = median(diff_boot);

    % Bootstrap percentile interval of observed p4-p3 differences
    diff_boot_CI = prctile(diff_boot, [100*alpha/2, 100*(1-alpha/2)]);

    % Bootstrap sign-based two-tailed p value
    % This asks whether the observed bootstrap difference distribution
    % is consistently above or below zero.
    n_left  = sum(diff_boot <= 0);
    n_right = sum(diff_boot >= 0);

    p_boot = 2 * min( ...
        (n_left  + 1) / (B + 1), ...
        (n_right + 1) / (B + 1) ...
    );
    p_boot = min(p_boot, 1);

    %% -----------------------------
    %  Paired label-swap permutation test
    %  -----------------------------
    % Under the null hypothesis, p3 and p4 labels are exchangeable
    % within each paired bootstrap repeat.
    null_mean = nan(nperm, 1);
    null_median = nan(nperm, 1);

    for iperm = 1:nperm

        swap_flag = rand(B, 1) > 0.5;

        perm_p1 = p1_boot;
        perm_p2 = p2_boot;

        tmp = perm_p1(swap_flag);
        perm_p1(swap_flag) = perm_p2(swap_flag);
        perm_p2(swap_flag) = tmp;

        perm_diff = perm_p2 - perm_p1;

        null_mean(iperm) = mean(perm_diff);
        null_median(iperm) = median(perm_diff);

    end

    % Two-tailed permutation p values
    p_perm_mean = (sum(abs(null_mean) >= abs(diff_mean)) + 1) / (nperm + 1);
    p_perm_median = (sum(abs(null_median) >= abs(diff_median)) + 1) / (nperm + 1);

    % Permutation null intervals
    % These are NOT confidence intervals of the observed difference.
    % They describe the expected null range after paired label swapping.
    null_mean_interval = prctile(null_mean, [100*alpha/2, 100*(1-alpha/2)]);
    null_median_interval = prctile(null_median, [100*alpha/2, 100*(1-alpha/2)]);

    %% -----------------------------
    %  Plot figure
    %  -----------------------------
    figure(pair_p(1)*100 + pair_p(2)*10 + iarea); clf;

    starttime_p1_val = cell2mat(starttime_p1(valid_p1));
    starttime_p2_val = cell2mat(starttime_p2(valid_p2));

    all_data = [starttime_p1_val(:); starttime_p2_val(:)];
    group = [ones(numel(starttime_p1_val),1); 2*ones(numel(starttime_p2_val),1)];

    hold on;
    colors = [p3_color; p4_color];

    boxplot(all_data, group, ...
        'orientation', 'horizontal', ...
        'Colors', colors, ...
        'Widths', 0.5, ...
        'Symbol', '', ...
        'Whisker', Inf);

    % Add scatter jitter along y-axis
    rng(1);
    scatter(starttime_p1_val(:), 1 + randn(numel(starttime_p1_val),1)*0.05, ...
        25, p3_color, 'filled');
    scatter(starttime_p2_val(:), 2 + randn(numel(starttime_p2_val),1)*0.05, ...
        25, p4_color, 'filled');

    yticks([1 2]);
    yticklabels({['p', num2str(pair_p(1))], ['p', num2str(pair_p(2))]});
    xlabel('Start time (ms)');
    xlim(xlim_range);
    ylim([0.5 2.5]);
    box off;
    set(gca, 'TickDir', 'out');

    title_info = sprintf('%s: p%d vs p%d, permutation p_{median}=%.5f', ...
        brain_area, pair_p(1), pair_p(2), p_perm_median);
    title(title_info, 'Interpreter', 'none');

    if exist('SetFigure', 'file') == 2
        SetFigure;
    end

    str = sprintf('p%d median = %.0f ms', pair_p(1), p1_median);
    annotation('textbox', [0.15 0.55 0.35 0.058], ...
        'String', str, 'EdgeColor', 'none', 'Color', p3_color);

    str = sprintf('p%d median = %.0f ms', pair_p(2), p2_median);
    annotation('textbox', [0.15 0.85 0.35 0.058], ...
        'String', str, 'EdgeColor', 'none', 'Color', p4_color);

    str = sprintf('median diff p%d-p%d = %.1f ms\nperm p = %.5f', ...
        pair_p(2), pair_p(1), diff_median, p_perm_median);
    annotation('textbox', [0.58 0.78 0.35 0.10], ...
        'String', str, 'EdgeColor', 'none');

    %% -----------------------------
    %  Print results
    %  -----------------------------
    fprintf('\n===== %s p%d start time =====\n', brain_area, pair_p(1));
    fprintf('Mean   = %.2f ms\n', p1_mean);
    fprintf('Median = %.2f ms\n', p1_median);
    fprintf('Bootstrap percentile interval = [%.2f, %.2f] ms\n', ...
        p1_boot_CI(1), p1_boot_CI(2));

    fprintf('\n===== %s p%d start time =====\n', brain_area, pair_p(2));
    fprintf('Mean   = %.2f ms\n', p2_mean);
    fprintf('Median = %.2f ms\n', p2_median);
    fprintf('Bootstrap percentile interval = [%.2f, %.2f] ms\n', ...
        p2_boot_CI(1), p2_boot_CI(2));

    fprintf('\n===== Bootstrap difference: p%d - p%d within %s =====\n', pair_p(2), pair_p(1), brain_area);
    fprintf('Mean difference across bootstrap samples   = %.2f ms\n', diff_mean);
    fprintf('Median difference across bootstrap samples = %.2f ms\n', diff_median);
    fprintf('Bootstrap 95%% interval of p%d-p%d difference = [%.2f, %.2f] ms\n', ...
        pair_p(2), pair_p(1), diff_boot_CI(1), diff_boot_CI(2));
    fprintf('Bootstrap sign-based p = %.5f\n', p_boot);

    fprintf('\n===== Paired label-swap permutation test =====\n');
    fprintf('Observed mean difference   = %.2f ms\n', diff_mean);
    fprintf('Permutation null 95%% interval for mean difference = [%.2f, %.2f] ms\n', ...
        null_mean_interval(1), null_mean_interval(2));
    fprintf('Permutation p for mean difference = %.5f\n', p_perm_mean);

    fprintf('\nObserved median difference = %.2f ms\n', diff_median);
    fprintf('Permutation null 95%% interval for median difference = [%.2f, %.2f] ms\n', ...
        null_median_interval(1), null_median_interval(2));
    fprintf('Permutation p for median difference = %.5f\n', p_perm_median);

    %% -----------------------------
    %  Optional interpretation
    %  -----------------------------
    fprintf('\n===== Interpretation guide =====\n');

    if diff_boot_CI(1) > 0
        fprintf('Bootstrap CI is entirely above 0: p%d starts later than p%d in %s.\n', ...
            pair_p(2), pair_p(1), brain_area);
    elseif diff_boot_CI(2) < 0
        fprintf('Bootstrap CI is entirely below 0: p%d starts earlier than p%d in %s.\n', ...
            pair_p(2), pair_p(1), brain_area);
    else
        fprintf('Bootstrap CI crosses 0: onset difference is not stable under bootstrap uncertainty.\n');
    end

    if p_perm_mean < alpha
        fprintf('Permutation test on mean difference is significant at alpha = %.2f.\n', alpha);
    else
        fprintf('Permutation test on mean difference is not significant at alpha = %.2f.\n', alpha);
    end

    if p_perm_median < alpha
        fprintf('Permutation test on median difference is significant at alpha = %.2f.\n', alpha);
    else
        fprintf('Permutation test on median difference is not significant at alpha = %.2f.\n', alpha);
    end

    %% -----------------------------
    %  Save results for this brain area
    %  -----------------------------
    R = struct();
    R.monkey = monkey;
    R.brain_area = brain_area;
    R.pair_p = pair_p;
    R.B = B;
    R.nperm = nperm;
    R.alpha = alpha;
    R.p1_boot = p1_boot;
    R.p2_boot = p2_boot;
    R.diff_boot = diff_boot;
    R.p1_mean = p1_mean;
    R.p2_mean = p2_mean;
    R.p1_median = p1_median;
    R.p2_median = p2_median;
    R.p1_boot_CI = p1_boot_CI;
    R.p2_boot_CI = p2_boot_CI;
    R.diff_mean = diff_mean;
    R.diff_median = diff_median;
    R.diff_boot_CI = diff_boot_CI;
    R.p_boot = p_boot;
    R.null_mean = null_mean;
    R.null_median = null_median;
    R.null_mean_interval = null_mean_interval;
    R.null_median_interval = null_median_interval;
    R.p_perm_mean = p_perm_mean;
    R.p_perm_median = p_perm_median;
    R.file_p1 = file_p1;
    R.file_p2 = file_p2;

    all_results.(brain_area) = R;

    if save_result
        save_file = ['./timedata/', monkey, brain_area, '_p', num2str(pair_p(1)), '_vs_p', num2str(pair_p(2)), '_Starttime_compare.mat'];
        save(save_file, 'R');
        fprintf('\nResults saved to:\n%s\n', save_file);
    end

end

if save_result
    save_file_all = ['./timedata/', monkey, '_p', num2str(pair_p(1)), '_vs_p', num2str(pair_p(2)), '_same_area_Starttime_compare_all.mat'];
    save(save_file_all, 'all_results');
    fprintf('\nAll-area results saved to:\n%s\n', save_file_all);
end

%% ========================================================================
%  Local function
%  ========================================================================
function starttime = local_load_starttime(file_name, brain_area)

    if ~exist(file_name, 'file')
        error('Input file does not exist: %s', file_name);
    end

    S = load(file_name);

    expected_var = [brain_area, '_starttime'];

    if isfield(S, expected_var)
        starttime = S.(expected_var);
    elseif isfield(S, 'starttime')
        starttime = S.starttime;
    else
        fields = fieldnames(S);
        idx = find(contains(lower(fields), 'starttime'), 1, 'first');
        if isempty(idx)
            error('Cannot find starttime variable in file: %s', file_name);
        else
            warning('Expected variable %s was not found. Using variable %s instead.', ...
                expected_var, fields{idx});
            starttime = S.(fields{idx});
        end
    end

end
