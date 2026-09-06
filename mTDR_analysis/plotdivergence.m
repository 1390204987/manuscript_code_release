%% 2025/2/9 calculate variable divergence
%% 2025/2/8 
% have 3 dfferent PSTH input, 
%1 cell input ,  cellsize = neuro_num*coherence_num*condition_num, under each cell matrixsizse = time_num*tiral_num
%2 matrix input, matrix size = neuro_num*time_num*coherence_num*condition_num
%3 matrix input, matrix size =
%neuro_num*time_num*trial_num*coherence_num*condition_num , for single trial analysis under simultaneous recording
%return diverge_period1PC1_conds, matrix size = time* reps 
function diverge_period1PC1_conds = plotdivergence(time,psth,unique_cons,control_cons,unique_cohe,Vec1,varargin)
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

        % 处理选择的 trials
        selecttrial_psth = cellfun(@(x) processCell(x, rep_num, time_num), ...
            psth, ...
            'UniformOutput', false);
        
        mean_psth = cellfun(@(x) nanmean(x, 2), psth, 'UniformOutput', false);
        
        
        for icohe = 1:length(unique_cohe)
            if isempty(control_cons)
                icon1 = find(unique_cons==-1);
                icon2 = find(unique_cons==1);
                temp1 = mean_psth(:,icohe,icon1)';
                temp2 = mean_psth(:,icohe,icon2)';
                mean_psth1 = cell2mat(temp1);
                mean_psth2 = cell2mat(temp2);
                mean_divergence = mean_psth1-mean_psth2;
                mean_divergence_period1PC1 = Vec1'*mean_divergence';
                
                psth1 = selecttrial_psth(:,icohe,icon1);
                psth2 = selecttrial_psth(:,icohe,icon2);
                reshaped_psth1 = cell2mat(cellfun(@(x) reshape(x, 1, []), psth1, 'UniformOutput', false));
                reshaped_psth2 = cell2mat(cellfun(@(x) reshape(x, 1, []), psth2, 'UniformOutput', false));
                selected_diverge = reshaped_psth1-reshaped_psth2 ;
                divergence_period1PC1 = Vec1'*selected_diverge;

                diverge_period1PC1 = reshape(divergence_period1PC1,[time_num,rep_num]); % 合并时间和条件维度          

            else
                diverge_period1PC1 = [];mean_divergence_period1PC1 = [];
                for icontrol_cons = 1:length(control_cons)% heading direction
                    icon1 = find(unique_cons==-1);
                    icon2 = find(unique_cons==1);
                    temp1 = mean_psth(:,icohe,icontrol_cons,icon1)';
                    temp2 = mean_psth(:,icohe,icontrol_cons,icon2)';
                    mean_psth1 = cell2mat(temp1);
                    mean_psth2 = cell2mat(temp2);
                    mean_divergence = mean_psth1-mean_psth2;
                    % for heading 0 codition
                    if length(unique_cons)==3
                        if icohe==1
                            heading_0 = find(unique_cons==0);
                            temp = mean_psth(:,icohe,icontrol_cons,heading_0)';
                            mean_psth0 = cell2mat(temp);
                            mean_divergence = mean_psth0;
                        end
                    end
                    temp_mean_diver_period1PC1{icontrol_cons} = Vec1'*mean_divergence';
                    
                    psth1 = selecttrial_psth(:,icohe,icontrol_cons,icon1);
                    psth2 = selecttrial_psth(:,icohe,icontrol_cons,icon2);
                    reshaped_psth1 = cell2mat(cellfun(@(x) reshape(x, 1, []), psth1, 'UniformOutput', false));
                    reshaped_psth2 = cell2mat(cellfun(@(x) reshape(x, 1, []), psth2, 'UniformOutput', false));
                    selected_diverge = reshaped_psth1-reshaped_psth2 ;
                    divergence_period1PC1= Vec1'*selected_diverge;
                    % for heading 0 codition
                    if length(unique_cons)==3
                        if icohe==1
                            heading_0 = find(unique_cons==0);
                            psth = selecttrial_psth(:,icohe,icontrol_cons,heading_0);
                            check = squeeze(selecttrial_psth(:,1,1,2));
                            reshaped_psth = cell2mat(cellfun(@(x) reshape(x, 1, []), psth, 'UniformOutput', false));
                            selected_diverge = reshaped_psth;
                            divergence_period1PC1 = Vec1'*selected_diverge;
                        end
                    end
                    temp_diver_period1PC1{icontrol_cons} = reshape(divergence_period1PC1,[time_num,rep_num]); % 合并时间和条件维度

                    if ~isnan( temp_mean_diver_period1PC1{icontrol_cons})
                        diverge_period1PC1 = [temp_diver_period1PC1{icontrol_cons},diverge_period1PC1];
                        mean_divergence_period1PC1 = [mean_divergence_period1PC1; temp_mean_diver_period1PC1{icontrol_cons}];
                    end
                end
                if size(mean_divergence_period1PC1,1)>1
                    mean_divergence_period1PC1 = mean(mean_divergence_period1PC1);
                end
            end
            this_color = min([color(icon1,:)+[1,1,1]*0.75*(1-1*icohe/length(unique_cohe));1,1,1]);
            all_diverge_period1PC1{icohe,1} = diverge_period1PC1;
            if icohe == length(unique_cohe)
                % 假设 all_diverge_period1PC1 是 3x1 的 cell，每个 cell 是 [T x N] 矩阵
                T = size(all_diverge_period1PC1{1}, 1);
                G = numel(all_diverge_period1PC1);
                
                % 构建一个函数，用于对每个时间点执行一次 ANOVA
                anova_p = @(t) ...
                    anova1( ...
                    cell2mat(cellfun(@(c) c(t, :)', all_diverge_period1PC1, 'UniformOutput', false)), ...
                    cell2mat(arrayfun(@(g) ...
                    repmat(g, size(all_diverge_period1PC1{g}(t, :), 2), 1), ...
                    (1:G)', 'UniformOutput', false)), ...
                    'off');
                
                % 用 arrayfun 执行每个时间点的分析
                p_values = arrayfun(anova_p, 1:T);
                sig_p_loc = p_values<0.001;
            end
            ys = mean_divergence_period1PC1;
            % resign the sign to prefer direction
            abs_value = abs(ys);
            maxindex = find(abs_value == max(abs_value));
            if ys(maxindex)>0
                ys = ys;
            else
                ys=-ys;
            end
            sem = nanstd(diverge_period1PC1, 0, 2)./sqrt(rep_num);
            CI = 1.96*sem;
            if ~isnan(ys)
                shadedErrorBar(time,ys,CI,'lineprops',{'Color',this_color,'LineStyle',lineStyle},'transparent',transparent);
                hold on
                if icohe == length(unique_cohe)
                    plot(time(sig_p_loc),zeros(size(time(sig_p_loc))),'*','Color',this_color)
                end
            end
%             plot(time,diverge_period1PC1(:,:,icohe),'Color',this_color,'LineWidth',linesize)
%             hold on
        end
        ylims = ylim;
        %     plot([seqtime1,seqtime1],ylims)
        xticks([0,1500]); xticklabels({'stimulus on ','stimulus off'})  
        

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

