%% 2025/2/9 calculate variable divergence
%% 2025/2/13
% have 3 dfferent PSTH input, 
%1 cell input ,  cellsize = neuro_num*coherence_num*condition_num, under each cell matrixsizse = time_num*tiral_num
%  or cellsize = neuro_num*coherence_num*condition1_num*condition2_num, under each cell matrixsizse = time_num*tiral_num
%2 matrix input, matrix size = neuro_num*time_num*coherence_num*condition_num
%3 matrix input, matrix size =
%neuro_num*time_num*trial_num*coherence_num*condition_num , for single trial analysis under simultaneous recording
%return diverge_period1PC1_conds, matrix size = time* reps 
function [mean_divergence,divergence] = getdivergence_new(time,psth,unique_cons,control_cons,unique_cohe,Vec1,coherence_0,varargin)
    if ~isempty(varargin)
        color = varargin{1};
        linesize = varargin{2};
        markersize = 10;
        transparent = 1;
        lineStyle = '-';
    else
        color = [1,0,0; 0,1,0];
        linesize = 2;
        markersize = 10;
        transparent = 1;
        lineStyle = '-';
    end
    time_num = size(time,2) ;

    if iscell(psth)
        if isempty(control_cons)
            rep_num = 50;% the number of selected trial each time
        else
            rep_num = 50;% the number of selected trial each time
        end
        if isempty(control_cons)
            [neuro_num,coh_num,con_num] = size(psth);
        else
            [neuro_num,coh_num,con1_num,con2_num] = size(psth);
            % con1 is heading con2 is choice
        end
      
        if coherence_0
            if isempty(control_cons)
                % 情况1：无 control_cons，仅处理 0 coherence (index=1)
                % 创建 neuron 和 condition 的索引 cell 数组

                [n_idx, s_idx] = ndgrid(1:neuro_num, 1:con_num);
                n_cell = num2cell(n_idx(:));
                s_cell = num2cell(s_idx(:));

                % 使用 cellfun 处理每个组合
                psth_conds = cellfun(@(n, s) horzcat(psth{n, 1, s}), ...
                    n_cell, s_cell, ...
                    'UniformOutput', false);
                % 处理选择的 trials
                select_conds_psth = cellfun(@(x) processCell(x, rep_num, time_num), ...
                    psth_conds, ...
                    'UniformOutput', false);
                
                % 如果需要保持原始维度 (neuro_num × con_num)
                psth_conds = reshape(psth_conds, neuro_num, con_num);
                select_conds_psth = reshape(select_conds_psth, neuro_num, con_num);

                select_cond_psth_mat = reshape(cell2mat(select_conds_psth),[time_num,neuro_num,rep_num,con_num]);
                select_cond_psth_mat  = permute(select_cond_psth_mat,[2,1,4,3]);       

                mean_psth = cellfun(@(x) nanmean(x, 2), psth_conds, 'UniformOutput', false);
            else
                [n_idx, s1_idx, s2_idx] = ndgrid(1:neuro_num, 1:con1_num,1:con2_num);
                n_cell = num2cell(n_idx(:));
                s1_cell = num2cell(s1_idx(:));
                s2_cell = num2cell(s2_idx(:));
                % 使用 cellfun 处理每个组合
                psth_conds = cellfun(@(n, s1, s2) horzcat(psth{n, 1, s1, s2}), ...
                                   n_cell, s1_cell, s2_cell, ...
                                   'UniformOutput', false);
               % 处理选择的 trials
               select_conds_psth = cellfun(@(x) processCell(x, rep_num, time_num), ...
                   psth_conds, ...
                   'UniformOutput', false);
                % 如果需要保持原始维度 (neuro_num ×coh_num x con1_num × con2_num)
                psth_conds = reshape(psth_conds, neuro_num,1, con1_num, con2_num);
                select_conds_psth = reshape(select_conds_psth, neuro_num, 1,con1_num, con2_num);
                mean_psth = cellfun(@(x) nanmean(x, 2), psth_conds, 'UniformOutput', false);
            end
            
        else
            if isempty(control_cons)
                % 创建神经元和条件的索引 cell 数组
                [n_idx, s_idx] = ndgrid(1:neuro_num, 1:con_num);
                n_cell = num2cell(n_idx(:));
                s_cell = num2cell(s_idx(:));
%                 % 使用 cellfun 合并所有 coherence 的数据
%                 psth_conds = cellfun(@(n, s) horzcat(psth{n, :, s}), ...
%                     n_cell, s_cell, ...
%                     'UniformOutput', false);
                
                % 处理选择的 trials
                select_conds_psth = cellfun(@(x) processCell(x, rep_num, time_num), ...
                    psth, ...
                    'UniformOutput', false);
                % 恢复原始维度 (neuro_num × con_num)
%                 psth_conds = reshape(psth_conds, neuro_num, con_num);

                select_conds_psth = reshape(select_conds_psth, neuro_num, coh_num, con_num);
%                 select_cond_psth_mat = reshape(cell2mat(select_conds_psth),[time_num,neuro_num,rep_num,con_num]);
%                 select_cond_psth_mat  = permute(select_cond_psth_mat,[2,1,4,3]);
                % here specificaly calculate mean_psth so we use all trial of
                % each neuro, or calculate from select_cond_psth_mat with
                % selectd trials
                mean_psth = cellfun(@(x) nanmean(x, 2), psth, 'UniformOutput', false);
            else               
                % 处理选择的 trials
                select_conds_psth = cellfun(@(x) processCell(x, rep_num, time_num), ...
                    psth, ...
                    'UniformOutput', false);
                
                % 恢复原始维度 (neuro_num x coh_num × con1_num × con2_num)

                select_conds_psth = reshape(select_conds_psth, neuro_num, coh_num,con1_num, con2_num);
                
                % 计算平均 PSTH
                mean_psth = cellfun(@(x) nanmean(x, 2), psth, 'UniformOutput', false);
                
            end
        end

        
        ylims = ylim;
        %     plot([seqtime1,seqtime1],ylims)
        xticks([0,1500]); xticklabels({'stimulus on ','stimulus off'})  
        if isempty(control_cons)
            mean_divergence = [];
            divergence = [];
            if coherence_0
                unique_cohe=1;
            end
            for icoh = 1:length(unique_cohe)
                icon1 = find(unique_cons==-1);
                icon2 = find(unique_cons==1);
                psth1 = select_conds_psth(:,icoh,icon1);
                psth2 = select_conds_psth(:,icoh,icon2);
                %             diverge_conds= squeeze(psth1_conds-psth2_conds);
                temp1 = mean_psth(:,icoh,icon1)';
                temp2 = mean_psth(:,icoh,icon2)';
                mean_psth1_conds = cell2mat(temp1);
                mean_psth2_conds = cell2mat(temp2);
                mean_diverge_conds = mean_psth1_conds-mean_psth2_conds;
                %             reshaped_diverge = reshape(diverge_conds, [neuro_num, time_num * rep_num]);
                %             rows_with_nan = any(isnan(reshaped_diverge),2);
                %             notnan_neuron = find(rows_with_nan==0);
                %             diverge_period1PC1_conds = Vec1'*reshaped_diverge(notnan_neuron,:)./sqrt(neuro_num);
                %             divergence = reshape(diverge_period1PC1_conds,[time_num, rep_num]);
                
                mat_psth1 = cat(3,psth1{:});
                mat_psth1 = permute(mat_psth1,[3,1,2]);
                mat_psth2 = cat(3,psth2{:});
                mat_psth2 = permute(mat_psth2,[3,1,2]);
                %too noisy for single trial divergence so here used randomly
                %pick and average method
                %设置参数
                num_repeats = 80;  % 总共挑选10次
                num_select = 30;    % 每次挑选5个重复试次（列）
%                 num_select = 50;
%                 % 使用 arrayfun 返回 cell 数组
%                 select_mat_psth1  = arrayfun(@(x) mat_psth1(:,:, randperm(rep_num, num_select)), ...
%                     1:num_repeats, 'UniformOutput', false);
%                 select_mat_psth2= arrayfun(@(x) mat_psth2(:,:, randperm(rep_num, num_select)), ...
%                     1:num_repeats, 'UniformOutput', false);
                num_neuron = size(mat_psth1, 1);
                num_time   = size(mat_psth1, 2);

                select_mat_psth1 = cell(1,num_repeats);
                select_mat_psth2 = cell(1,num_repeats);
                base_seed = 20260427;  % 你自己指定一个固定种子
                for r = 1:num_repeats
                    tmp1 = nan(num_neuron, num_time, num_select);
                    tmp2 = nan(num_neuron, num_time, num_select);
                    for n = 1:num_neuron
                        rng(base_seed + r+n, 'twister');  % 每一次 repeat 使用特定随机数种子
                        perm1 = randperm(rep_num, num_select);
                        perm2 = randperm(rep_num, num_select);
                        tmp1(n,:,:) = mat_psth1(n,:,perm1);
                        tmp2(n,:,:) = mat_psth2(n,:,perm2);
                    end
                    select_mat_psth1{r} = tmp1;
                    select_mat_psth2{r} = tmp2;
                end
                
                reshaped_psth1 = cellfun(@(x) reshape(x, neuro_num, []), select_mat_psth1, 'UniformOutput', false);
                reshaped_psth2 = cellfun(@(x) reshape(x, neuro_num, []), select_mat_psth2, 'UniformOutput', false);
                psth1_period1PC1 = cellfun(@(x) Vec1'*x./sqrt(neuro_num),reshaped_psth1,'UniformOutput', false) ;
                psth2_period1PC1 = cellfun(@(x) Vec1'*x./sqrt(neuro_num),reshaped_psth2,'UniformOutput', false) ;
                psth1_period1PC1 = cellfun(@(x) reshape(x,[time_num,num_select]),psth1_period1PC1,'UniformOutput', false) ;
                psth2_period1PC1 = cellfun(@(x) reshape(x,[time_num,num_select]),psth2_period1PC1,'UniformOutput', false) ;
                temp_mean_divergence{icoh} =  Vec1'*mean_diverge_conds'./sqrt(neuro_num);
                temp_divergence{icoh}  = cellfun(@(x,y) x-y, psth1_period1PC1,psth2_period1PC1,'UniformOutput', false);
                
                if ~isnan(temp_mean_divergence{icoh})
                    divergence = [temp_divergence{icoh},divergence];
                    mean_divergence = [mean_divergence;temp_mean_divergence{icoh}];
                end
            end

        else
            mean_divergence = [];
            divergence = [];
            if coherence_0
                unique_cohe=1;
            end
            for icoh = 1:length(unique_cohe)
                for icontrol_cons = 1:length(control_cons) % heading direction
                    icon1 = find(unique_cons==-1);
                    icon2 = find(unique_cons==1);
                    psth1 = select_conds_psth(:,icoh,icontrol_cons,icon1);
                    psth2 = select_conds_psth(:,icoh,icontrol_cons,icon2);
                    mean_psth1 =  mean_psth(:,icoh,icontrol_cons,icon1);
                    mean_psth2 = mean_psth(:,icoh,icontrol_cons,icon2);
                    mean_psth1_mat = cell2mat(mean_psth1');
                    mean_psth2_mat = cell2mat(mean_psth2');
                    mat_psth1 = cat(3,psth1{:});
                    mat_psth1 = permute(mat_psth1,[3,1,2]);
                    mat_psth2 = cat(3,psth2{:});
                    mat_psth2 = permute(mat_psth2,[3,1,2]);
                    %too noisy for single trial divergence so here used randomly
                    %pick and average method
                    %设置参数
                    num_repeats = 50;  % 总共挑选10次
                    num_select = 30;    % 每次挑选5个重复试次（列）
%                     num_select = 50;  
%                     % 使用 arrayfun 返回 cell 数组
%                     select_mat_psth1  = arrayfun(@(x) mat_psth1(:,:, randperm(rep_num, num_select)), ...
%                         1:num_repeats, 'UniformOutput', false);
                    num_neuron = size(mat_psth1, 1);
                    num_time   = size(mat_psth1, 2);

                    select_mat_psth1 = cell(1,num_repeats);
                    select_mat_psth2 = cell(1,num_repeats);
                    base_seed = 20260427;
                    for r = 1:num_repeats
                        tmp1 = nan(num_neuron, num_time, num_select);
                        tmp2 = nan(num_neuron, num_time, num_select);
                        for n = 1:num_neuron
                            rng(base_seed + r+n, 'twister');
                            perm1 = randperm(rep_num, num_select);
                            perm2 = randperm(rep_num, num_select);
                            tmp1(n,:,:) = mat_psth1(n,:,perm1);
                            tmp2(n,:,:) = mat_psth2(n,:,perm2);
                        end
                        select_mat_psth1{r} = tmp1;
                        select_mat_psth2{r} = tmp2;
                    end
%                     select_mat_psth2= arrayfun(@(x) mat_psth2(:,:, randperm(rep_num, num_select)), ...
%                         1:num_repeats, 'UniformOutput', false);
%                     psth1_period1PC1 = cellfun(@(x) mean(x, 2), selected_psth1_period1PC1, 'UniformOutput', false);
%                     psth2_period1PC1 = cellfun(@(x) mean(x, 2), selected_psth2_period1PC1, 'UniformOutput', false);
%                     temp_divergence{icontrol_cons,icoh} = cell2mat(psth1_period1PC1)-cell2mat(psth2_period1PC1);

                    reshaped_psth1 = cellfun(@(x) reshape(x, neuro_num, []), select_mat_psth1, 'UniformOutput', false);
                    reshaped_psth2 = cellfun(@(x) reshape(x, neuro_num, []), select_mat_psth2, 'UniformOutput', false);
                                       
%                     if size(Vec1,2)>1
%                         % get main axis of pojection  data pca
%                         temp1 = Vec1'*reshaped_psth1;
%                         temp2 = Vec1'*reshaped_psth2;
%                         [PCs_psth1,~,~] = pca(temp1');
%                         [PCs_psth2,~,~] = pca(temp2');
%                         pc3_psth1 = PCs_psth1(:,:);
%                         pc3_psth2 = PCs_psth2(:,:);
%                         % denoise of projection data
%                         temp1 =pc3_psth1*pc3_psth1'*temp1;
%                         temp2 =pc3_psth2*pc3_psth2'*temp2;
%                         psth1_period1PC1 = sum(temp1)./sqrt(neuro_num);
%                         psth2_period1PC1 = sum(temp2)./sqrt(neuro_num);
%                         % same process on mean psth
%                         temp1 = Vec1'*mean_psth1_mat';
%                         temp2 = Vec1'*mean_psth2_mat';
%                         [PCs_psth1,~,~] = pca(temp1');
%                         [PCs_psth2,~,~] = pca(temp2');
%                         pc3_psth1 = PCs_psth1(:,:);
%                         pc3_psth2 = PCs_psth2(:,:);
%                         temp1 =pc3_psth1*pc3_psth1'*temp1;
%                         temp2 =pc3_psth2*pc3_psth2'*temp2;
%                         mean_psth1_period1PC1 = sum(temp1)./sqrt(neuro_num);
%                         mean_psth2_period1PC1 = sum(temp2)./sqrt(neuro_num);
%                     else
                        psth1_period1PC1 = cellfun(@(x) Vec1'*x./sqrt(neuro_num),reshaped_psth1,'UniformOutput', false) ;
                        psth2_period1PC1 = cellfun(@(x) Vec1'*x./sqrt(neuro_num),reshaped_psth2,'UniformOutput', false) ;
                        mean_psth1_period1PC1 = Vec1'*mean_psth1_mat'./sqrt(neuro_num);
                        mean_psth2_period1PC1 = Vec1'*mean_psth2_mat'./sqrt(neuro_num);
%                     end
                    psth1_period1PC1 = cellfun(@(x) reshape(x,[time_num,num_select]),psth1_period1PC1,'UniformOutput', false) ;
                    psth2_period1PC1 = cellfun(@(x) reshape(x,[time_num,num_select]),psth2_period1PC1,'UniformOutput', false) ;
%                     temp_divergence{icontrol_cons,icoh} = psth1_period1PC1-psth2_period1PC1;
                    temp_mean_divergence{icontrol_cons,icoh} = mean_psth1_period1PC1-mean_psth2_period1PC1;                    
                    temp_divergence{icontrol_cons,icoh}  = cellfun(@(x,y) x-y, psth1_period1PC1,psth2_period1PC1,'UniformOutput', false);
                    if ~isnan(temp_mean_divergence{icontrol_cons,icoh})
                        divergence = [temp_divergence{icontrol_cons,icoh},divergence];
                        mean_divergence = [mean_divergence;temp_mean_divergence{icontrol_cons,icoh}];
                    end
                end
            end
        end
        if coherence_0==0
            mean_divergence = mean(mean_divergence);
        end
    else
        for icohe = 1:length(unique_cohe)
            icon1 = find(unique_cons==-1);
            icon2 = find(unique_cons==1);
            psth1 = psth(:,:,icohe,icon1);
            psth2 = psth(:,:,icohe,icon2);
            diverge(:,:,icohe) = squeeze(psth1-psth2);
            diverge_period1PC1(:,:,icohe) = Vec1'*diverge(:,:,icohe);
            this_color = min([color(icon1,:)+[1,1,1]*0.75*(1-1*icohe/length(unique_cohe));1,1,1]);
            plot(time,diverge_period1PC1(:,:,icohe),'Color',this_color,'LineWidth',linesize)
            hold on
        end
        ylims = ylim;
        %     plot([seqtime1,seqtime1],ylims)
        xticks([0,1500]); xticklabels({'stimulus on ','stimulus off'})  
        mean_diverge_period1PC1 = nanmean(diverge_period1PC1,3);
        std_diverge_period1PC1 = nanstd(diverge_period1PC1,0,3);       
        this_color = min([color(icon1,:)+[1,1,1]*0.75*(1-1*1);1,1,1]);
        ys= mean_diverge_period1PC1;
        sem=std_diverge_period1PC1/sqrt(length(unique_cohe));
        CI = sem*1.96;
        if ~isnan(ys)
            shadedErrorBar(time,ys,CI,'lineprops',{'Color',this_color,'LineStyle',lineStyle},'transparent',transparent);
            hold on
        end
        diverge_period1PC1 = squeeze(diverge_period1PC1);
    end
    
    function processedCell = processCell(current_psth,rep_num,time_num)
        if ~isempty(current_psth)
            current_reps = size(current_psth, 2); % 当前 trial 数量
            
            if current_reps < rep_num
                % 如果 trial 数量不足 rep_num，随机重复 trial
                idx = randi(current_reps, 1, rep_num); % 随机选择 trial 索引（允许重复）
                processedCell = current_psth(:, idx); % 抽取或重复 trial
%                 % 保留已有 trial，不足的补 NaN
%                 processedCell = nan(time_num, rep_num);
%                 processedCell(:, 1:current_reps) = current_psth;
            elseif current_reps > rep_num
                % 如果 trial 数量超过 rep_num，随机抽取 5 个 trial
                idx = randperm(current_reps, rep_num); % 随机选择 5 个 trial 索引
                processedCell = current_psth(:, idx); % 抽取 trial
            else
                % 如果 trial 数量正好是 5，直接使用
                processedCell = current_psth;
            end
        else
            processedCell  = nan(time_num,rep_num);% 39is the time points
        end
    end
        
    end

