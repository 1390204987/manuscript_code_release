% this is for get shuffle actvitiy divergence , skip random selection to
% same time
function [mean_divergence] = getdivergence_shuffle(time,psth,unique_cons,control_cons,unique_cohe,Vec1,coherence_0,varargin)
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

                
                % 如果需要保持原始维度 (neuro_num × con_num)
                psth_conds = reshape(psth_conds, neuro_num, con_num);
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
                % 如果需要保持原始维度 (neuro_num ×coh_num x con1_num × con2_num)
                psth_conds = reshape(psth_conds, neuro_num,1, con1_num, con2_num);
                mean_psth = cellfun(@(x) nanmean(x, 2), psth_conds, 'UniformOutput', false);
            end
        else
            if isempty(control_cons)
                mean_psth = cellfun(@(x) nanmean(x, 2), psth, 'UniformOutput', false);
            else     
                % 计算平均 PSTH
                mean_psth = cellfun(@(x) nanmean(x, 2), psth, 'UniformOutput', false);
            end
        end
        ylims = ylim;
        %     plot([seqtime1,seqtime1],ylims)
        xticks([0,1500]); xticklabels({'stimulus on ','stimulus off'})  
        if isempty(control_cons)
            mean_divergence = [];
            if coherence_0
                unique_cohe=1;
            end
            for icoh = 1:length(unique_cohe)
                icon1 = find(unique_cons==-1);
                icon2 = find(unique_cons==1);

                %             diverge_conds= squeeze(psth1_conds-psth2_conds);
                temp1 = mean_psth(:,icoh,icon1)';
                temp2 = mean_psth(:,icoh,icon2)';
                mean_psth1_conds = cell2mat(temp1);
                mean_psth2_conds = cell2mat(temp2);
                mean_diverge_conds = mean_psth1_conds-mean_psth2_conds;
                
                temp_mean_divergence{icoh} =  Vec1'*mean_diverge_conds'./sqrt(neuro_num);
                
                if ~isnan(temp_mean_divergence{icoh})
                    mean_divergence = [mean_divergence;temp_mean_divergence{icoh}];
                end
            end
        else
            mean_divergence = [];
            
            if coherence_0
                unique_cohe=1;
            end
            for icoh = 1:length(unique_cohe)
                for icontrol_cons = 1:length(control_cons) % heading direction
                    icon1 = find(unique_cons==-1);
                    icon2 = find(unique_cons==1);
                    mean_psth1 =  mean_psth(:,icoh,icontrol_cons,icon1);
                    mean_psth2 = mean_psth(:,icoh,icontrol_cons,icon2);
                    mean_psth1_mat = cell2mat(mean_psth1');
                    mean_psth2_mat = cell2mat(mean_psth2');
                    mean_psth1_period1PC1 = Vec1'*mean_psth1_mat'./sqrt(neuro_num);
                    mean_psth2_period1PC1 = Vec1'*mean_psth2_mat'./sqrt(neuro_num);
                    temp_mean_divergence{icontrol_cons,icoh} = mean_psth1_period1PC1-mean_psth2_period1PC1;                    
                    if ~isnan(temp_mean_divergence{icontrol_cons,icoh})
                        mean_divergence = [mean_divergence;temp_mean_divergence{icontrol_cons,icoh}];
                    end
                end
            end
        end
        if coherence_0==0
            mean_divergence = mean(mean_divergence);
        end
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