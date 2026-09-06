%% compare_starttime_bootstrap_integrated.m
% Directly load saved MST/LIP start-time files.
% Compare onset/start time between MST and LIP using:
% 1) Bootstrap difference distribution
% 2) Paired label-swap permutation test
%
% Important:
% - Bootstrap CI is computed from observed bootstrap differences:
%       diff_boot = LIP_boot - MST_boot
% - Permutation null interval is computed from shuffled paired labels.
%   It is a null interval, not the confidence interval of the observed difference.

clear; clc;

%% -----------------------------
%  Parameters
%  -----------------------------
monkey = 'F';
p = 3;                  % modify this
nperm = 5000;           % permutation repeats
alpha = 0.05;

MST_file = ['./timedata/',monkey,'p', num2str(p), 'MSTtime.mat'];
LIP_file = ['./timedata/',monkey,'p', num2str(p), 'LIPtime.mat'];

save_file = ['./timedata/p', num2str(p), 'Starttime_compare_integrated.mat'];
color_{1,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
color_{2,1} = [0,191,255;155,48,255]/255; %color
color_{3,1} = [139,10,80;72,61,139]/255; %choice
color_{4,1} = [0,0,128;205,133,0]/255; %saccade
color_{5,1} = [0,0,139;205,133,0]/255; %last saccade
color_{6,1} = [0,0,1;0,0,0.9;0,0,0.8;0,0,0.7;0,0,0.6]; %coherence
color_{7,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
color = color_{p,1};
% this_color = min([color(1,:)+[1,1,1]*1*(1-1*ibrain_area/4);1,1,1]);

if p==4
    this_color = [30,114,255]/255;
    MST_color = this_color;
else
    MST_color = min([color(1,:)+[1,1,1]*1*(1-1*2/4);1,1,1]);
end


if p==4
    this_color = [30,64,139]/255;
    LIP_color = this_color;
else
    LIP_color = min([color(1,:)+[1,1,1]*1*(1-1*4/4);1,1,1]);
end



%% -----------------------------
%  Load saved start-time files
%  -----------------------------
if ~exist(MST_file, 'file')
    error('MST file does not exist: %s', MST_file);
end

if ~exist(LIP_file, 'file')
    error('LIP file does not exist: %s', LIP_file);
end

load(MST_file, "MST_starttime");
load(LIP_file, "LIP_starttime");

%% -----------------------------
%  Make sure data are cell arrays
%  -----------------------------
if ~iscell(MST_starttime)
    MST_starttime = num2cell(MST_starttime);
end

if ~iscell(LIP_starttime)
    LIP_starttime = num2cell(LIP_starttime);
end

%% -----------------------------
%  Remove empty bootstrap samples
%  -----------------------------
valid_MST = ~cellfun(@isempty, MST_starttime);
valid_LIP = ~cellfun(@isempty, LIP_starttime);

% Paired bootstrap comparison:
% only keep repeats where both MST and LIP have valid onset times.
valid_pair = valid_MST & valid_LIP;

MST_boot = cell2mat(MST_starttime(valid_pair));
LIP_boot = cell2mat(LIP_starttime(valid_pair));

MST_boot = MST_boot(:);
LIP_boot = LIP_boot(:);

B = length(MST_boot);

if B < 5
    error('Too few valid paired bootstrap samples: B = %d', B);
end

if length(MST_boot) ~= length(LIP_boot)
    error('MST_boot and LIP_boot have different lengths.');
end

fprintf('\nNumber of valid paired bootstrap samples = %d\n', B);

%% -----------------------------
%  Descriptive statistics for each area
%  -----------------------------
MST_mean = mean(MST_boot);
LIP_mean = mean(LIP_boot);

MST_median = median(MST_boot);
LIP_median = median(LIP_boot);

% Since MST_boot and LIP_boot are already bootstrap distributions,
% these percentile ranges describe bootstrap uncertainty of each area's onset.
MST_boot_CI = prctile(MST_boot, [100*alpha/2, 100*(1-alpha/2)]);
LIP_boot_CI = prctile(LIP_boot, [100*alpha/2, 100*(1-alpha/2)]);

%% -----------------------------
%  Bootstrap difference distribution
%  -----------------------------
% Positive value:
%   LIP starts later than MST
%
% Negative value:
%   LIP starts earlier than MST

diff_boot = LIP_boot - MST_boot;

diff_mean = mean(diff_boot);
diff_median = median(diff_boot);

% Bootstrap percentile interval of observed LIP-MST differences
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
% Under the null hypothesis, MST and LIP labels are exchangeable
% within each paired bootstrap repeat.

null_mean = nan(nperm, 1);
null_median = nan(nperm, 1);

for iperm = 1:nperm

    swap_flag = rand(B, 1) > 0.5;

    perm_MST = MST_boot;
    perm_LIP = LIP_boot;

    tmp = perm_MST(swap_flag);
    perm_MST(swap_flag) = perm_LIP(swap_flag);
    perm_LIP(swap_flag) = tmp;

    perm_diff = perm_LIP - perm_MST;

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

figure(p*10+2)
MST_starttime_val = cell2mat(MST_starttime);
LIP_starttime_val = cell2mat(LIP_starttime);
all_data = [cell2mat(MST_starttime), cell2mat(LIP_starttime)];
group = [ones(size(MST_starttime_val)), 2*ones(size(LIP_starttime_val))];
hold on
colors = [MST_color;LIP_color];
boxplot(all_data', group','orientation','horizontal', 'Colors', colors, 'Widths', 0.5, 'Symbol', '','Whisker', Inf);
% 添加散点 jitter（横向 jitter 是 y轴扰动）
rng(1);
scatter(cell2mat(MST_starttime), 1 + randn(numel(MST_starttime_val),1)*0.05, 25, MST_color, 'filled');  % 蓝色
scatter(cell2mat(LIP_starttime), 2 + randn(numel(LIP_starttime_val),1)*0.05, 25, LIP_color, 'filled');    % 橙色               
% 美化
yticks([1 2]); yticklabels({'Group 1', 'Group 2'});
xlabel('Response');
%                 xlim([0,1000])
xlim([200,1300])
ylim([0.5 2.5]);
box off; set(gca, 'TickDir', 'out');
title_info = ['p=',num2str(p_perm_median)];
title(title_info)
SetFigure
%                 str = sprintf('MST median = %d\n' , median(cell2mat(MST_starttime)));
%                 annotation('textbox', [0.15 0.85 0.35 0.058], 'String', str);
%                 str = sprintf('LIP median = %d\n' , median(cell2mat(LIP_starttime)));
%                 annotation('textbox', [0.15 0.55 0.35 0.058], 'String', str);
str = sprintf('MST median = %d\n' , MST_median);
annotation('textbox', [0.15 0.55 0.35 0.058], 'String', str);
str = sprintf('LIP median = %d\n' , LIP_median );
annotation('textbox', [0.15 0.85 0.35 0.058], 'String', str);
%% -----------------------------
%  Print results
%  -----------------------------
fprintf('\n===== MST start time =====\n');
fprintf('Mean   = %.2f ms\n', MST_mean);
fprintf('Median = %.2f ms\n', MST_median);
fprintf('Bootstrap percentile interval = [%.2f, %.2f] ms\n', ...
    MST_boot_CI(1), MST_boot_CI(2));

fprintf('\n===== LIP start time =====\n');
fprintf('Mean   = %.2f ms\n', LIP_mean);
fprintf('Median = %.2f ms\n', LIP_median);
fprintf('Bootstrap percentile interval = [%.2f, %.2f] ms\n', ...
    LIP_boot_CI(1), LIP_boot_CI(2));

fprintf('\n===== Bootstrap difference: LIP - MST =====\n');
fprintf('Mean difference across bootstrap samples   = %.2f ms\n', diff_mean);
fprintf('Median difference across bootstrap samples = %.2f ms\n', diff_median);
fprintf('Bootstrap 95%% interval of LIP-MST difference = [%.2f, %.2f] ms\n', ...
    diff_boot_CI(1), diff_boot_CI(2));
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
    fprintf('Bootstrap CI is entirely above 0: LIP starts later than MST under bootstrap uncertainty.\n');
elseif diff_boot_CI(2) < 0
    fprintf('Bootstrap CI is entirely below 0: LIP starts earlier than MST under bootstrap uncertainty.\n');
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

% %% -----------------------------
% %  Save results
% %  -----------------------------
% save(save_file, ...
%     "p", "B", "nperm", "alpha", ...
%     "MST_boot", "LIP_boot", "diff_boot", ...
%     "MST_mean", "LIP_mean", ...
%     "MST_median", "LIP_median", ...
%     "MST_boot_CI", "LIP_boot_CI", ...
%     "diff_mean", "diff_median", ...
%     "diff_boot_CI", "p_boot", ...
%     "null_mean", "null_median", ...
%     "null_mean_interval", "null_median_interval", ...
%     "p_perm_mean", "p_perm_median");
% 
% fprintf('\nResults saved to:\n%s\n', save_file);