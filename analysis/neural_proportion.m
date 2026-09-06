%calculate neural proportion encoding saccade choice , perceptual choice or
%both
close all; clear;
%%
% FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseSTS_LA\'};
% FolderPath_HD =  {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorSTS_LA\'; 'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorSTS_NP\'};
FolderPath_HD =  {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorNTIPS_LA\'; 'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTIPS_NP\'};
% 
% FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseIPS_LA\'};
% FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorSTS_LA\'; 'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorSTS_NP\';...
%     'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseSTS_LA\'};
if length(FolderPath_HD)==1
    Files_HD = dir(fullfile(FolderPath_HD{1},'*PSTH.mat'));
    pathfilenum = length(Files_HD);
else
    Files_HD = [];
    for ipath = 1:length(FolderPath_HD)
        ipathfiles = dir(fullfile(FolderPath_HD{ipath},'*PSTH.mat'));
        Files_HD = [Files_HD;ipathfiles];
        pathfilenum(ipath) = length(ipathfiles);
    end
    pathfilenum = cumsum(pathfilenum);
end
File_work_count = 0;
for ct = 1:length(Files_HD)
    File_work_count = File_work_count+1;
    current_path_pod = find(pathfilenum<File_work_count);
    if isempty(current_path_pod)
        current_path = 1;
    else
        current_path = max(current_path_pod) + 1;
    end
    list{File_work_count} = load([FolderPath_HD{current_path} Files_HD(ct).name]);
    cell_channel = list{1,ct}.result.SpikeChan;
    TF = isstrprop(list{1,ct}.result.FILE,'digit');
    charloc = find(TF==0);
    monkeyID(ct) = str2num(list{1,ct}.result.FILE(charloc(1)+1:charloc(2)-1));
    cellnum(ct) = str2double(list{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1));
    cellID(ct) = monkeyID(ct)*10^6+str2double(list{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^3+cell_channel;
end

%%  select specific cells(have significant 1D tuning or mem sac tuning)

%cellID_1D
%     FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dMST_LA\','Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dMT_LA\'};
% FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\Friday1dSTS_NP\'};
FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dIPS_LA\'};
% FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\Friday1dSTS_NP\'};
%     FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dIPS_LA\'};
%     FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_NP\Friday1dIPS_NP\'};

% FolderPath_1D = {'Z:\Data\Tempo\Batch\test\Friday1d112_NP\'};
% FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_NP\Friday1dNTSTS_NP\'};
if length(FolderPath_1D)==1
    Files_1D = dir(fullfile(FolderPath_1D{1},'*AzimuthTuning.mat'));
    pathfilenum_1D = length(Files_1D);
else
    Files_1D = [];
    for ipath = 1:length(FolderPath_1D)
        ipathfiles = dir(fullfile(FolderPath_1D{ipath},'*AzimuthTuning.mat'));
        Files_1D = [Files_1D;ipathfiles];
        pathfilenum_1D(ipath) = length(ipathfiles);
    end
    pathfilenum_1D = cumsum(pathfilenum_1D);
end
File_work_count_1D = 0;
for ct = 1:length(Files_1D)
    File_work_count_1D = File_work_count_1D+1;
    current_path_pod = find(pathfilenum_1D<File_work_count_1D);
    if isempty(current_path_pod)
        current_path = 1;
    else
        current_path = max(current_path_pod) + 1;
    end
    list_1D{File_work_count_1D} = load([FolderPath_1D{current_path} Files_1D(ct).name]);
    cell_channel_1D = list_1D{1,ct}.result.SpikeChan;
    p_1D(ct) = list_1D{1,ct}.result.p_1D;
    az_1D(ct) = list_1D{1,ct}.result.az;
    TF = isstrprop(list_1D{1,ct}.result.FILE,'digit');
    charloc = find(TF==0);
    monkeyID_1D = str2num(list_1D{1,ct}.result.FILE(charloc(1)+1:charloc(2)-1));
    cellnum_1D(ct) = str2double(list_1D{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1));
    cellID_1D(ct) = str2double(list_1D{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel_1D;
end
select_cellID_1D = cellID_1D(p_1D<1);

%%
%     %cellID_mem
%    FolderPath_mem = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogonMemIPS_LA\'};
FolderPath_mem = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaymemSTS_LA\'};
if length(FolderPath_mem)==1
    Files_mem = dir(fullfile(FolderPath_mem{1},'*MemSac.mat'));
    pathfilenum = length(Files_mem);
else
    Files_mem = [];
    pathfilenum = [];
    for ipath = 1:length(FolderPath_mem)
        ipathfiles = dir(fullfile(FolderPath_mem{ipath},'*MemSac.mat'));
        Files_mem = [Files_mem;ipathfiles];
        pathfilenum(ipath) = length(ipathfiles);
    end
    pathfilenum = cumsum(pathfilenum);
end

File_work_count_mem = 0;
for ct = 1:length(Files_mem)
    File_work_count_mem = File_work_count_mem+1;
    current_path_pod = find(pathfilenum<File_work_count_mem);
    if isempty(current_path_pod)
        current_path = 1;
    else
        current_path = max(current_path_pod) + 1;
    end
    list_mem{File_work_count_mem} = load([FolderPath_mem{current_path} Files_mem(ct).name]);
    p_mem(ct) = min(list_mem{1,ct}.result.p);
    cell_channel = list_mem{1,ct}.result.SpikeChan;
    TF = isstrprop(list_mem{1,ct}.result.FILE,'digit');
    charloc = find(TF==0);
    monkeyID = str2num(list_mem{1,ct}.result.FILE(charloc(1)+1:charloc(2)-1));
    cellnum_mem(ct) = str2double(list_mem{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;
    blocknum_mem(ct) = str2double(list_mem{1,ct}.result.FILE(charloc(3)+1:end));
    cellID_mem(ct) = str2double(list_mem{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;   %+blocknum(ct)/10;
end
select_cellID_mem = cellID_mem(p_mem<0.05);
%%
% [cellID_pick,cellindex,cellindex_1D] = intersect(cellID,select_cellID_1D);
%     [cellID_pick,cellindex,cellindex_mem] = intersect(cellID,select_cellID_mem);
[cellID_pick,cellindex,cellindex] = intersect(cellID,cellID);

neuralnum = length(cellID_pick)
%%  
% set time related parameter
%time window
p_critical = 0.05;
align_offsets_others = list{1,1}.result.align_offsets_others;
stim_dur = align_offsets_others{1,2}(1,2) - align_offsets_others{1,2}(1,1);
stimulus_on = 0;
stimulus_end = 0+stim_dur;
%method1 calcalate with psth data
ts_per_trial = list{1,1}.result.t_centers{1,1}' ;
time_window1(1,1) = find(ts_per_trial(:) > stimulus_on-800,1);
time_window1(1,2) = find(ts_per_trial(:) < stimulus_end,1,'last');
time_window1(1,3) = size(ts_per_trial,1);
%method2 calculate with 0/1 spike data
spike_aligned = list{1,1}.result.spike_aligned;
time_window2(1,1) = find(spike_aligned{2,1}(:) > stimulus_on+500,1);
time_window2(1,2) = find(spike_aligned{2,1}(:) < stimulus_end,1,'last');
time_window2(1,3) = size(spike_aligned{2,1},2);
%%
for f = 1:length(cellID_pick)
    heading_per_trial = (list{1,cellindex(f)}.result.heading_per_trial)';
    coherence_per_trial = (list{1,cellindex(f)}.result.coherence_per_trial)';
    T1loc_per_trial = (list{1,cellindex(f)}.result.targetloc_per_trial)';
    choice_per_trial = (list{1,cellindex(f)}.result.choice_per_trial)';
    sac_per_trial = sign(T1loc_per_trial).*choice_per_trial; %
%     spike_hist = list{1,cellindex(f)}.result.spike_hist{1,1};
    spike_aligned = list{1,cellindex(f)}.result.spike_aligned;
    align_markers = (list{1,cellindex(f)}.result.align_markers);
    align_offsets_others = list{1,cellindex(f)}.result.align_offsets_others;

    if length(unique(coherence_per_trial))==1
        relativeheading_per_trial = heading_per_trial;
    else
        relativeheading_per_trial = sign(heading_per_trial).*coherence_per_trial;
    end
    unirelativeheading = unique(relativeheading_per_trial);
    unichoice = unique(choice_per_trial);
% calculate firing rate
    for j =1
        fir_per_trial = sum(spike_aligned{1,j}(:,time_window2(1,1):time_window2(1,2)),2)/((time_window2(1,2)-time_window2(1,1))/1000);
    end
   % roc calculate perceptual choice selectivity
   for iheading = 1:length(unirelativeheading)
       select = (relativeheading_per_trial == unirelativeheading(iheading));
       z_heading =  fir_per_trial(select,:);
       z_heading_temp = (z_heading - mean(z_heading,1))/std(z_heading(:));
       Zspike_heading_pertrial(select,:) = z_heading_temp;
   end
   Zresp_left_choice =  Zspike_heading_pertrial(choice_per_trial==-1,:);
   Zresp_right_choice = Zspike_heading_pertrial(choice_per_trial==1,:);
   rng(1);  % repeatable random Nan20260408
   permuteN=1000;
    [~,shuffle_index] = sort(rand(length(Zspike_heading_pertrial),permuteN),1); % This is really brilliant.
   if length(Zresp_left_choice) >3 && length(Zresp_right_choice)>3
       [choice_roc(f,1),~,perm] = rocN(Zresp_left_choice,Zresp_right_choice,100,1000);  %set the left as x axis
      for ishuffle = 1:1000
          Zspike_heading_pertrial_shuff = Zspike_heading_pertrial(shuffle_index(:,ishuffle),:);
          Zresp_left_choice_shuff =  Zspike_heading_pertrial_shuff(choice_per_trial ==-1,:);
          Zresp_right_choice_shuff = Zspike_heading_pertrial_shuff(choice_per_trial ==1,:);
          [shuffle_choice(f,ishuffle),~,~] = rocN(Zresp_left_choice_shuff,Zresp_right_choice_shuff,100);
      end
%       pvalue_choice_roc(f,1) = sum(shuffle_choice(f,:)>choice_roc(f,1))/1000;
%       if choice_roc(f,1)<0.5, pvalue_choice_roc(f,1)  = 1-pvalue_choice_roc(f,1); end % One tail
%       pvalue_choice_roc(f,1) = pvalue_choice_roc(f,1)*2; % Two tail
          pvalue_choice_roc(f,1) = perm.pValue ;
   else
       [choice_roc(f,1),~,perm] = NaN;
   end
   % roc calculate saccade choice selectivity
   Zresp_left_sac =  Zspike_heading_pertrial(sac_per_trial ==-1,:);
   Zresp_right_sac = Zspike_heading_pertrial(sac_per_trial ==1,:);
      if length(Zresp_left_sac) >3 && length(Zresp_right_sac)>3
          [sac_roc(f,1),~,perm] = rocN(Zresp_left_sac,Zresp_right_sac,100,1000);  %set the left as x axis
          for ishuffle = 1:1000
              Zspike_heading_pertrial_shuff = Zspike_heading_pertrial(shuffle_index(:,ishuffle),:);
              Zresp_left_sac_shuff =  Zspike_heading_pertrial_shuff(sac_per_trial ==-1,:);
              Zresp_right_sac_shuff = Zspike_heading_pertrial_shuff(sac_per_trial ==1,:);
              [shuffle_sac(f,ishuffle),~,~] = rocN(Zresp_left_sac_shuff,Zresp_right_sac_shuff,100);
          end
%           pvalue_sac_roc(f,1) = sum(shuffle_sac(f,:)>sac_roc(f,1))/1000;
%           if sac_roc(f,1)<0.5, pvalue_sac_roc(f,1)  = 1-pvalue_sac_roc(f,1); end % One tail
%           pvalue_sac_roc(f,1) = pvalue_sac_roc(f,1)*2; % Two tail
              pvalue_sac_roc(f,1)  = perm.pValue;
      end
%    if length(Zresp_left_sac) >3 && length(Zresp_right_sac)>3
%        [sac_roc(f,1),~,perm] = rocN(Zresp_left_sac,Zresp_right_sac,100,1000);  %set the left as x axis
%        pvalue_sac_roc(f,1) = perm.pValue; 
%        if ~isnan(sac_roc(f,1))
%            shuffle_sac(f,:) = perm.auROCPerm';
%        else
%            shuffle_sac(f,:)  = nan(1,1000);
%        end
% %        [~,pvalue_sac_roc(f,1)] = ttest2(Zresp_left_sac,Zresp_right_sac);
%    else
%        [sac_roc(f,1),~,perm] = NaN;
%    end
end
% 需要下载venn函数（可从MATLAB File Exchange获取）
% 示例：https://www.mathworks.com/matlabcentral/fileexchange/27076-venn
sig_abschoice = sum(pvalue_choice_roc<0.05&pvalue_sac_roc>=0.05); 
sig_sac = sum(pvalue_sac_roc<0.05&pvalue_choice_roc>=0.05);
sig_both = sum(pvalue_choice_roc<0.05&pvalue_sac_roc<0.05);
check = sum(pvalue_choice_roc>=0.05&pvalue_sac_roc>=0.05)
%% statistical check
%{
%-----------------------binomial test-------------------------------------------------
% 1. 参数设置
N = length(pvalue_choice_roc); % 总细胞数
alpha = 0.05;                  % 单个神经元的显著性水平
% 理论上的 Chance Level 概率
p_chance_single = alpha;               % 对于单一维度 (Choice 或 Saccade)
p_chance_both   = alpha * alpha;       % 对于交集 (假设两者独立，0.0025)
% 2. 进行二项分布检验 (Binomial Test)
% binocdf(k, N, p) 计算的是观察到 <=k 次成功的概率
% 我们需要计算观察到 >=k 次成功的概率，即 1 - binocdf(k-1, N, p)
% --- Choice 显著数量检验 ---
p_val_abschoice = 1 - binocdf((sig_abschoice+sig_both) - 1, N, p_chance_single);
% --- Saccade 显著数量检验 ---
p_val_sac = 1 - binocdf((sig_sac+sig_both) - 1, N, p_chance_single);
% --- Both (两者皆显著) 数量检验 ---
p_sac_empirical = sig_sac / N; % 使用实测的 Saccade 假阳性率
p_val_both = 1 - binocdf(sig_both - 1, sig_abschoice+sig_both,p_sac_empirical);
% 3. 结果汇总与显示
stats_labels = {'Choice ', 'Saccade', 'Both Significant'};
stats_counts = [(sig_abschoice+sig_both), (sig_sac+sig_both), sig_both];
stats_pvals  = [p_val_abschoice, p_val_sac, p_val_both];
stats_chances = [N * p_chance_single, N * p_chance_single, N * p_chance_both];

fprintf('\n--- 显著性细胞统计与 Chance Level 检验 ---\n');
fprintf('%-20s | %-10s | %-15s | %-10s\n', 'Category', 'Count', 'Chance Level', 'P-value');
fprintf('----------------------------------------------------------------------\n');

for i = 1:3
    res_str = ' ';
    if stats_pvals(i) < 0.05, res_str = '* (Sig)'; end
    if stats_pvals(i) < 0.01, res_str = '** (Highly Sig)'; end
    
    fprintf('%-20s | %-10d | %-15.2f | %.4f %s\n', ...
        stats_labels{i}, stats_counts(i), stats_chances(i), stats_pvals(i), res_str);
end
%}
%--------------------------------------permute test--------------------------------
n_neurons = size(shuffle_choice, 1);
n_shuffles = size(shuffle_choice, 2);
alpha = 0.05;                  % 单个神经元的显著性水平
% 理论上的 Chance Level 概率
p_chance_single = alpha;               % 对于单一维度 (Choice 或 Saccade)
p_chance_both   = alpha * alpha;       % 对于交集 (假设两者独立，0.0025)
% --- 第一步：构建每个神经元的显著性判断矩阵 ---
sig_mat_choice = false(n_neurons, n_shuffles);
sig_mat_sac = false(n_neurons, n_shuffles);

for n = 1:n_neurons
    % 为每个神经元单独计算阈值 (双尾)
    th_choice = prctile(shuffle_choice(n,:), [2.5, 97.5]);
    sig_mat_choice(n,:) = shuffle_choice(n,:) < th_choice(1) | shuffle_choice(n,:) > th_choice(2);
    th_sac = prctile(shuffle_sac(n,:), [2.5, 97.5]);
    sig_mat_sac(n,:) = shuffle_sac(n,:) < th_sac(1) | shuffle_sac(n,:) > th_sac(2);
end

% --- 第二步：计算 Choice-Only 和 Sac-Only 的随机分布 ---
% Choice-only 的定义：Choice 显著 且 Saccade 不显著
null_choice_counts = sum(sig_mat_choice, 1); 
% Saccade-only 的定义：Saccade 显著 且 Choice 不显著
null_sac_counts = sum(sig_mat_sac, 1);
% null both 的定义： choice显著且 saccade 显著
null_both = sum(sig_mat_choice.*sig_mat_sac, 1);

% --- 第三步：统计各成分细胞数量的显著性 (P-value) ---
pvalue_choice_N= sum(null_choice_counts >= (sig_abschoice+sig_both)) / n_shuffles;
pvalue_sac_N = sum(null_sac_counts >= (sig_sac+sig_both)) / n_shuffles;

pvalue_onlychoice_N= sum(null_choice_counts >= (sig_abschoice)) / n_shuffles;
pvalue_onlysac_N = sum(null_sac_counts >= (sig_sac)) / n_shuffles;

% --- 第四步：统计interleave细胞数量的显著性 (P-value) ---
p_choice = (pvalue_choice_roc < 0.05); % 逻辑向量
p_sac = (pvalue_sac_roc < 0.05);
real_both = sum(p_choice & p_sac);
% p_perm_fast = sum(null_both >= real_both) / 1000;
% chance_both = prctile(null_both,95);

n_shufs = 1000;
null_both = zeros(n_shufs, 1);
for i = 1:n_shufs
    null_both(i) = sum(p_choice & p_sac(randperm(length(p_sac))));
end
p_perm_fast = sum(null_both >= real_both) / n_shufs;
chance_both = prctile(null_both,95);


% --- 第五步：统计两种细胞数量差异的显著性 (P-value) ---
diff_null_sac_choice = -(null_choice_counts - null_sac_counts);
real_diff_sac_choice = -(sig_abschoice - sig_sac);
% p_diff_sac_choice = mean(diff_null_sac_choice >= real_diff_sac_choice);
% chance_diff_sac_choice = prctile(abs(diff_null_sac_choice),95);
% two tailed
p_diff_sac_choice = mean(abs(diff_null_sac_choice) >= abs(real_diff_sac_choice));
chance_diff_sac_choice = prctile(abs(diff_null_sac_choice),95);

diff_null_sac_both = null_sac_counts - null_both;
real_diff_sac_both = sig_sac - sig_both;
% p_diff_sac_both = mean(diff_null_sac_both>= real_diff_sac_both);
% chance_diff_sac_both = prctile(abs(diff_null_sac_both),95);
% two tailed
p_diff_sac_both = mean(abs(diff_null_sac_both) >= abs(real_diff_sac_both));
chance_diff_sac_both = prctile(abs(diff_null_sac_both), 95);

stats_labels = {'Choice ', 'Saccade', 'Both Significant', 'diff sac abs Sig', 'diff sac both Sig'};
stats_counts = [(sig_abschoice), (sig_sac), sig_both, real_diff_sac_choice, real_diff_sac_both];
stats_pvals  = [pvalue_choice_N, pvalue_sac_N, p_perm_fast, p_diff_sac_choice, p_diff_sac_both];
stats_chances = [n_neurons * p_chance_single, n_neurons* p_chance_single,chance_both,chance_diff_sac_choice, chance_diff_sac_both];

fprintf('\n--- 显著性细胞统计与 Chance Level 检验 ---\n');
fprintf('%-20s | %-10s | %-15s | %-10s\n', 'Category', 'Count', 'Chance Level', 'P-value');
fprintf('----------------------------------------------------------------------\n');

for i = 1:5
    res_str = ' ';
    if stats_pvals(i) < 0.05, res_str = '* (Sig)'; end
    if stats_pvals(i) < 0.01, res_str = '** (Highly Sig)'; end
    
    fprintf('%-20s | %-10d | %-15.2f | %.4f %s\n', ...
        stats_labels{i}, stats_counts(i), stats_chances(i), stats_pvals(i), res_str);
end

%%
% index_abs = pvalue_choice_roc<0.05;
% list_abs = list(1,index_abs);
% popu_CP(list_abs)
%%
% 创建精确数值的韦恩图
colorabs = [139,10,80]/255; %choice
colorsac = [0,0,139]/255; %saccade
alpha_value = 0.5;

% 计算圆半径（面积与数据大小成比例）
total = sig_abschoice + sig_sac + sig_both;
r_choice = sqrt(sig_abschoice + sig_both)/sqrt(total)*0.8;
r_sac = sqrt(sig_sac + sig_both)/sqrt(total)*0.8;

% 计算圆心位置（确保重叠区域面积正确）
d = sqrt(r_choice^2 + r_sac^2 - 2*sqrt(sig_both/total)*r_choice*r_sac);

% venn([sig_abschoice, sig_sac, sig_both], ...
%     'FaceColor', {colorabs, colorsac}, ...      % 红色代表choice，蓝色代表sac
%     'FaceAlpha', {0.5, 0.5}, ...     % 50%透明度
%     'EdgeColor', 'k', ...            % 黑色边框
%     'LineWidth', 2)           % 边框宽度);           
% % 直接在图中标注数值
% text(-0.5, 0, num2str(sig_abschoice), 'FontSize', 12, 'HorizontalAlignment', 'center');
% text(9, 0, num2str(sig_sac), 'FontSize', 12, 'HorizontalAlignment', 'center');
% text(5, 0, num2str(sig_both), 'FontSize', 12, 'HorizontalAlignment', 'center');

% 创建自定义韦恩图
figure;
hold on;

% 绘制choice圆
rectangle('Position', [-d/2-r_choice, -r_choice, 2*r_choice, 2*r_choice], ...
          'Curvature', [1 1], 'FaceColor', [colorabs,alpha_value], 'EdgeColor', 'k', 'LineWidth', 2);

% 绘制sac圆
rectangle('Position', [d/2-r_sac, -r_sac, 2*r_sac, 2*r_sac], ...
          'Curvature', [1 1], 'FaceColor', [colorsac,alpha_value], 'EdgeColor', 'k', 'LineWidth', 2);

% 添加数值标签
text(-d/2, 0, num2str(sig_abschoice), 'FontSize', 12, 'HorizontalAlignment', 'center');
text(d/2, 0, num2str(sig_sac), 'FontSize', 12, 'HorizontalAlignment', 'center');
text(0, 0, num2str(sig_both), 'FontSize', 12, 'HorizontalAlignment', 'center');

% 添加标题和图例
title(['Significant Signals Overlap (Choice: ', num2str(sig_abschoice), ...
       ' | SAC: ', num2str(sig_sac), ' | Both: ', num2str(sig_both), ')']);
% legend({'Choice only', 'SAC only', 'Both'}, 'Location', 'best');
xlim([-1.1,1.1])
ylim([-1.1,1.1])
% axis equal
% axis off
% SetFigure
% axis equal tight  % 同时保证坐标轴比例和紧凑显示
set(gcf, 'Position', [100 100 600 600])  % 设置方形图形窗口
axis off