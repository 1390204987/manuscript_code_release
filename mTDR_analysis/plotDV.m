%% 2025/2/23 calculate unsigned decision variable
%% 2025/2/8 
% have 3 dfferent PSTH input, 
%1 cell input ,  cellsize = neuro_num*coherence_num*condition_num, under each cell matrixsizse = time_num*tiral_num
%2 matrix input, matrix size = neuro_num*time_num*coherence_num*condition_num
%3 matrix input, matrix size =
%neuro_num*time_num*trial_num*coherence_num*condition_num , for single trial analysis under simultaneous recording
%return diverge_period1PC1_conds, matrix size = time* reps 
function DV_period1PC1_conds = plotDV(time,psth,unique_cons,unique_cohe,Vec1,varargin)
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
        rep_num = 15;% the number of selected trial each time
        [neuro_num,coh_num,con_num] = size(psth);


        selecttrial_psth = arrayfun(@(n, c, s) processCell(psth{n, c, s},rep_num,time_num), ...
            repmat((1:neuro_num)', 1, coh_num, con_num), ... % 神经元索引
            repmat(1:coh_num, neuro_num, 1, con_num), ...    % 一致性索引
            repmat(reshape(1:con_num, 1, 1, con_num), neuro_num, coh_num, 1), ... % 条件索引
            'UniformOutput', false);


        mean_psth = cellfun(@(x) nanmean(x, 2), psth, 'UniformOutput', false);
        
        % 将元胞数组转换为矩阵
        temp = cell2mat(selecttrial_psth); % 转换为 (neuro_num*time_num) x (coh_num*repetition_num) x con_num
        % 调整矩阵形状
        temp = reshape(temp, time_num, neuro_num, rep_num, coh_num, con_num); % 调整为 time_num x neuro_num x coh_num x con_num x repetition_num
        
        % 调整维度顺序
        selecttrial_psth_mat = permute(temp, [2, 1, 4, 5, 3]); % 调整为 neuro_num x time_num x coh_num x con_num x repetition_num
        

        for icohe = 1:length(unique_cohe)
            icon1 = find(unique_cons==-1);
            icon2 = find(unique_cons==1);
            psth1 = squeeze(selecttrial_psth_mat(:,:,icohe,icon1,:));
            psth2 = squeeze(selecttrial_psth_mat(:,:,icohe,icon2,:));
            temp1 = mean_psth(:,icohe,icon1)';
            temp2 = mean_psth(:,icohe,icon2)';
            mean_psth1 = cell2mat(temp1);
            mean_psth2 = cell2mat(temp2);
            % 提取数据并调整维度
            reshape_psth1 =  reshape(psth1, [neuro_num, time_num * rep_num]); % 合并时间和条件维度
            reshape_psth2 =  reshape(psth2, [neuro_num, time_num * rep_num]); % 合并时间和条件维度
            %             rows_with_nan = any(isnan(reshaped_diverge),2);
            %             notnan_neuron = find(rows_with_nan==0);
            % 执行矩阵乘法
            %             result = Vec1(notnan_neuron)' * reshaped_diverge(notnan_neuron,:); % 结果维度: 1 × (time_num * con_num)
            psth1_period1PC1 = Vec1'*reshape_psth1;
            psth2_period1PC1 = Vec1'*reshape_psth2;
            % 恢复原始时间和条件维度
            psth1_period1PC1  = reshape(psth1_period1PC1, [time_num, rep_num]); % 维度: time_num × con_num
            psth2_period1PC1  = reshape(psth2_period1PC1, [time_num, rep_num]); % 维度: time_num × con_num
            mean_psth1_period1PC1 = Vec1'*mean_psth1';
            mean_psth2_period1PC1 = Vec1'*mean_psth2';
            this_color = min([color(icon1,:)+[1,1,1]*0.75*(1-1*icohe/length(unique_cohe));1,1,1]);
            unsignedpsth_period1PC1 = [psth1_period1PC1,-psth2_period1PC1];
            mean_unsignedpsth_period1PC1 = mean([mean_psth1_period1PC1;-mean_psth2_period1PC1]);
%             unsignedpsth_period1PC1 = psth1_period1PC1;
%             ys = mean(unsignedpsth_period1PC1,2);
            ys = mean_unsignedpsth_period1PC1;
            % resign the sign to prefer direction
            abs_value = abs(ys);
            maxindex = find(abs_value == max(abs_value));
            if ys(maxindex)>0
                ys = ys;
            else
                ys=ys;
            end
            sem = nanstd(unsignedpsth_period1PC1, 0, 2)./sqrt(rep_num);
            CI = 1.96*sem;
            if ~isnan(ys)
%                 shadedErrorBar(time,ys,CI,'lineprops',{'Color',this_color,'LineStyle',lineStyle},'transparent',transparent);
                plot(time,ys,'Color',this_color,'LineWidth',linesize)
                hold on
            end
            %             plot(time,diverge_period1PC1(:,:,icohe),'Color',this_color,'LineWidth',linesize)
            %             hold on
        end
        ylims = ylim;
    end
    
    function processedCell = processCell(current_psth,rep_num,time_num)
        if ~isempty(current_psth)
            current_reps = size(current_psth, 2); % 当前 trial 数量
            
            if current_reps < rep_num
                % 如果 trial 数量不足 rep_num，随机重复 trial
                idx = randi(current_reps, 1, rep_num); % 随机选择 trial 索引（允许重复）
                processedCell = current_psth(:, idx); % 抽取或重复 trial
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

