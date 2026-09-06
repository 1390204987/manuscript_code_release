%% 2025/2/8 
% have 3 dfferent PSTH input, 
%1 cell input ,  cellsize = neuro_num*coherence_num*condition_num, under each cell matrixsizse = time_num*tiral_num
%2 matrix input, matrix size = neuro_num*time_num*coherence_num*condition_num
%3 matrix input, matrix size =
%neuro_num*time_num*trial_num*coherence_num*condition_num , for single trial analysis under simultaneous recording
function plotpsth1D(time,psth,unique_cons,unique_cohe,Vec1,varargin)
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
        rep_num = 15;
        [neuro_num,coh_num,con_num] = size(psth);
%         reps = cellfun(@(x) size(x,2), psth);
        % 使用 arrayfun 处理每个 cell
        selecttrial_psth = arrayfun(@(n, c, s) processCell(psth{n, c, s},rep_num,time_num), ...
            repmat((1:neuro_num)', 1, coh_num, con_num), ... % 神经元索引
            repmat(1:coh_num, neuro_num, 1, con_num), ...    % 一致性索引
            repmat(reshape(1:con_num, 1, 1, con_num), neuro_num, coh_num, 1), ... % 条件索引
            'UniformOutput', false);
        
        mean_psth = arrayfun(@(n, c, s) nanmean(selecttrial_psth{n, c, s}, 2), ...
            repmat((1:neuro_num)', 1, coh_num, con_num), ... % 神经元索引
            repmat(1:coh_num, neuro_num, 1, con_num), ...    % 一致性索引
            repmat(reshape(1:con_num, 1, 1, con_num), neuro_num, coh_num, 1), ... % 条件索引
            'UniformOutput', false);
        
        std_psth = arrayfun(@(n, c, s) nanstd(selecttrial_psth{n, c, s}, 0, 2), ...
            repmat((1:neuro_num)', 1, coh_num, con_num), ... % 神经元索引
            repmat(1:coh_num, neuro_num, 1, con_num), ...    % 一致性索引
            repmat(reshape(1:con_num, 1, 1, con_num), neuro_num, coh_num, 1), ... % 条件索引
            'UniformOutput', false);
        
        % 将 mean_psth 和 std_psth 转换为目标大小的矩阵
        mean_psth_matrix = reshape(cell2mat(mean_psth), [time_num, neuro_num, coh_num, con_num]);
        mean_psth_matrix = permute(mean_psth_matrix, [2, 1, 3, 4]); % 调整为 neuro_num x time_num x coh_num x con_num
        
        std_psth_matrix = reshape(cell2mat(std_psth), [time_num, neuro_num, coh_num, con_num]);
        std_psth_matrix = permute(std_psth_matrix, [2, 1, 3, 4]); % 调整为 neuro_num x time_num x coh_num x con_num
%         for icons = 1:length(unique_cons)
%             for icohe = 1:length(unique_cohe)
%                 mean_psth_select = mean_psth_matrix(:,:,icohe,icons);
%                 rows_with_nan = any(isnan(mean_psth_select),2);
%                 notnan_neuron = find(rows_with_nan==0);
%                 mean_psth_period1PC1 = Vec1(notnan_neuron)'*mean_psth_select(notnan_neuron,:);               
%                 ys = mean_psth_period1PC1;
%                 std_psth_select = std_psth_matrix(:,:,icohe,icons);
%                 std_psth_period1PC1 = Vec1(notnan_neuron)'*std_psth_select(notnan_neuron,:);
%                 sem = std_psth_period1PC1./sqrt(rep_num);
%                 CI = sem*1.96;            
%                 this_color = min([color(icons,:)+[1,1,1]*0.75*(1-1*icohe/length(unique_cohe));1,1,1]);
%                 index_start = find(time>-500,1);
%                 index_end = find(time<2000,1,'last');
%                 if ~isnan(ys)
%                     shadedErrorBar(time,ys,CI,'lineprops',{'Color',this_color,'LineStyle',lineStyle},'transparent',transparent);
%                     hold on
%                 end
%             end
%         end
        for icons = 1:length(unique_cons)
            for icohe = 1:length(unique_cohe)

                % -------------------------------------------------------------
                % Reconstruct trial-level population PSTH
                % pop_psth: neuron_num × time_num × rep_num
                % -------------------------------------------------------------
                pop_psth = nan(neuro_num, time_num, rep_num);

                for n = 1:neuro_num
                    this_cell = selecttrial_psth{n, icohe, icons};  % time_num × rep_num

                    if ~isempty(this_cell)
                        pop_psth(n,:,:) = permute(this_cell, [3, 1, 2]);
                    end
                end

                % -------------------------------------------------------------
                % Remove neurons with NaN values
                % -------------------------------------------------------------
                rows_with_nan = squeeze(any(any(isnan(pop_psth), 2), 3));
                notnan_neuron = find(rows_with_nan == 0);

                if isempty(notnan_neuron)
                    continue
                end

                % -------------------------------------------------------------
                % Trial-by-trial projection onto Vec1
                % projected_trial: rep_num × time_num
                % -------------------------------------------------------------
                projected_trial = nan(rep_num, time_num);

                for r = 1:rep_num
                    this_trial_psth = pop_psth(notnan_neuron,:,r);  % neuron × time

                    projected_trial(r,:) = Vec1(notnan_neuron)' * this_trial_psth;

                    % 如果你想除以 sqrt(neuron number)，用下面这行替代上一行：
                    % projected_trial(r,:) = Vec1(notnan_neuron)' * this_trial_psth ./ sqrt(length(notnan_neuron));
                end

                % -------------------------------------------------------------
                % Compute mean and error bar across projected trials
                % -------------------------------------------------------------
                ys = nanmean(projected_trial, 1);                       % 1 × time_num
                sem = nanstd(projected_trial, 0, 1) ./ sqrt(rep_num);    % 1 × time_num
                CI = 1.96 .* sem;                                       % 95% CI approximation

                % -------------------------------------------------------------
                % Plot
                % -------------------------------------------------------------
                this_color = min([
                    color(icons,:) + [1,1,1] * 0.75 * (1 - icohe / length(unique_cohe));
                    1,1,1
                ]);

                if any(~isnan(ys))
                    shadedErrorBar(time, ys, CI, ...
                        'lineprops', {'Color', this_color, 'LineStyle', lineStyle}, ...
                        'transparent', transparent);
                    hold on
                end

            end
        end

        ylims = ylim;
        %     plot([seqtime1,seqtime1],ylims)
        xticks([0,1500]); xticklabels({'stimulus on ','stimulus off'})
        
    else
                
        for icons = 1:length(unique_cons)
            for icohe = 1:length(unique_cohe)
                psth_select = psth(:,:,icohe,icons);
                psth_period1PC1 = Vec1'*psth_select;
                this_color = min([color(icons,:)+[1,1,1]*0.75*(1-1*icohe/length(unique_cohe));1,1,1]);
                index_start = find(time>-500,1);
                index_end = find(time<2000,1,'last');
                %             plot(time(index_start:index_end),psth_period1PC1(index_start:index_end),'Color',this_color,'LineWidth',linesize)
                plot(time,psth_period1PC1,'Color',this_color,'LineWidth',linesize)
                hold on
            end
        end
        ylims = ylim;
        %     plot([seqtime1,seqtime1],ylims)
        xticks([0,1500]); xticklabels({'stimulus on ','stimulus off'})
%         yticks(ylims)
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