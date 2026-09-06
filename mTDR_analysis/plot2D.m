%% 2025/2/9 
% have 3 dfferent PSTH input, 
%1 cell input ,  cellsize = neuro_num*coherence_num*condition_num, under each cell matrixsizse = time_num*tiral_num
%2 matrix input, matrix size = neuro_num*time_num*coherence_num*condition_num
%3 matrix input, matrix size =
%neuro_num*time_num*trial_num*coherence_num*condition_num , for single trial analysis under simultaneous recording

function plot2D(time,psth,unique_cons,unique_cohe,Vec1,Vec2,varargin)

    if ~isempty(varargin)
        color = varargin{1};
        linesize = varargin{2};
        markersize = 10;
    else
        color = [1,0,0; 0,1,0];
        linesize = 2;
        markersize = 10;
        
    end
    time_num = size(time,2) ;
    if iscell(psth)
        rep_num = 10;
        [neuro_num,coh_num,con_num] = size(psth);
        reps = cellfun(@(x) size(x,2), psth);
        reps = reshape(reps,[neuro_num,coh_num,con_num]);
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
        
        for icons = 1:length(unique_cons)
            for icohes = 1:length(unique_cohe)
                psth_select = mean_psth_matrix(:,:,icohes,icons);
                psth_Vec1 = Vec1'*psth_select;
                psth_Vec2 = Vec2'*psth_select;
                %             this_color = color(icons,:)*(1-0.3*icohes/length(unique_cohe))+[0,0,1]*0.3*icohes/length(unique_cohe);
                this_color = min([color(icons,:)+[1,1,1]*0.75*(1-1*icohes/length(unique_cohe));1,1,1]);
                index_start = find(time>-500,1);
                index_end = find(time<2500,1,'last');
                plot(psth_Vec1(index_start:index_end),psth_Vec2(index_start:index_end),'Color',this_color,'LineWidth',linesize)
                hold on
                index_time0 = find(time>0,1);
                index_stim_off = find(time<1500,1,'last');
                index_750 = find(time>750,1);
                plot(psth_Vec1(index_time0),psth_Vec2(index_time0),'o','Color',this_color,'MarkerFaceColor',this_color,'MarkerSize',markersize)
                plot(psth_Vec1(index_stim_off),psth_Vec2(index_stim_off),'^','Color',this_color,'MarkerFaceColor',this_color,'MarkerSize',markersize)
                plot(psth_Vec1(index_750),psth_Vec2(index_750),'o','Color',this_color,'MarkerFaceColor',this_color,'MarkerSize',5)
            end
        end

    else
        for icons = 1:length(unique_cons)
            for icohes = 1:length(unique_cohe)
                psth_select = psth(:,:,icohes,icons);
                psth_Vec1 = Vec1'*psth_select;
                psth_Vec2 = Vec2'*psth_select;
                %             this_color = color(icons,:)*(1-0.3*icohes/length(unique_cohe))+[0,0,1]*0.3*icohes/length(unique_cohe);
                this_color = min([color(icons,:)+[1,1,1]*0.75*(1-1*icohes/length(unique_cohe));1,1,1]);
                index_start = find(time>-500,1);
                index_end = find(time<2500,1,'last');
                plot(psth_Vec1(index_start:index_end),psth_Vec2(index_start:index_end),'Color',this_color,'LineWidth',linesize)
                hold on
                index_time0 = find(time>0,1);
                index_stim_off = find(time<1500,1,'last');
                index_750 = find(time>750,1);
                plot(psth_Vec1(index_time0),psth_Vec2(index_time0),'o','Color',this_color,'MarkerFaceColor',this_color,'MarkerSize',markersize)
                plot(psth_Vec1(index_stim_off),psth_Vec2(index_stim_off),'^','Color',this_color,'MarkerFaceColor',this_color,'MarkerSize',markersize)
                plot(psth_Vec1(index_750),psth_Vec2(index_750),'o','Color',this_color,'MarkerFaceColor',this_color,'MarkerSize',5)
            end
        end
        
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