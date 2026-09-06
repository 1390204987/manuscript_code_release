% changed the shuffle 
%% by zn 2025/2/11 compare variable encoding in different brain area
close all;clear;
% MT_data = load('D_MTcolor_all');
% MT_vec = load('D_colorvariable_MT.mat');

monkey = 'D';
% monkey = 'F';
MT_data = [];
MT_seqPCA_vec = [];
MT_mTDR_vec = [];

MST_data = load([monkey,'_STScolor_all_20ms']);
MST_seqPCA_vec = load([monkey,'_colorvariable_all_STS.mat']);
MST_mTDR_vec = load([monkey,'_STScolor_all_mTDR.mat']);

% shu_MST_data = [];
% shu_MST_vec = [];
% MST_data = [];
% MST_vec = [];
% VIP_data = load('D_VIPcolor_all');
% VIP_vec = load('D_colorvariable_VIP.mat');

VIP_data =[];
VIP_seqPCA_vec = [];
VIP_mTDR_vec = [];

LIP_data =load([monkey,'_IPScolor_all_20ms']);
LIP_seqPCA_vec = load([monkey,'_colorvariable_all_IPS.mat']);
LIP_mTDR_vec = load([monkey,'_IPScolor_all_mTDR.mat']);

% LIP_data =[];
% LIP_vec = [];
% shu_LIP_data =[];
% shu_LIP_vec = [];
color_{1,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
color_{2,1} = [0,191,255;155,48,255]/255; %color
color_{3,1} = [139,10,80;72,61,139]/255; %choice
color_{4,1} = [0,0,128;205,133,0]/255; %saccade
color_{5,1} = [0,0,139;205,133,0]/255; %last saccade
color_{6,1} = [0,0,1;0,0,0.9;0,0,0.8;0,0,0.7;0,0,0.6]; %coherence
color_{7,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
linesize = 2;
transparent = 1;
lineStyle = '-';
variable_name={'heading';'color';'choice';'saccade';'last_sac';'coherence'};
p_critical = 0.001;
for p = 3:4
    color = color_{p,1};
    
    MST_starttime = [];
    LIP_starttime = [];
    for ibrain_area = 1:4 %1=MT,2=MST,3=VIP,4=LIP
        if ibrain_area == 1
            data =MT_data; vec_seqPCA = MT_seqPCA_vec; vec_mTDR = MT_mTDR_vec; 
        elseif ibrain_area == 2
            data =MST_data; vec_seqPCA = MST_seqPCA_vec; vec_mTDR = MST_mTDR_vec;               
        elseif ibrain_area == 3
            data =VIP_data; vec_seqPCA = VIP_seqPCA_vec; vec_mTDR = VIP_mTDR_vec;
        elseif ibrain_area == 4
            data =LIP_data; vec_seqPCA = LIP_seqPCA_vec; vec_mTDR = LIP_mTDR_vec; 
        end  
        if ~isempty(data)
            time = data.t_centers{:};
            if p == 1
                [psth,heading_cons,choice_cons,unique_cohe] = get_2v_PSTH(data.Xmat,data.Ymat,data.hk);
                unique_cons = heading_cons;
                control_cons = choice_cons;
                psth = permute(psth,[1,2,4,3]);
            elseif p == 3
                [psth,heading_cons,choice_cons,unique_cohe] = get_2v_PSTH(data.Xmat,data.Ymat,data.hk);
                unique_cons = choice_cons;
                control_cons = heading_cons;
            else
                [psth,~,~,unique_cons,unique_cohe] = get_PSTH(data.Xmat,data.Ymat,p,data.hk);
                control_cons = [];
            end
            icon1 = find(unique_cons==-1);
            icon2 = find(unique_cons==1);
            coherence_0=0; 

%             [mean_divergepsth_seqPC1,divergepsth_seqPC1]=getdivergence(time,psth,unique_cons,control_cons,unique_cohe,vec_seqPCA.variablePC{1,p}(:,1),coherence_0,color,linesize);
            [mean_divergepsth_seqPC1,divergepsth_seqPC1]=getdivergence_new(time,psth,unique_cons,control_cons,unique_cohe,vec_seqPCA.variablePC{1,p}(:,1),coherence_0,color,linesize);
            if ibrain_area == 2
                brain_area = 'STS';
            elseif ibrain_area == 4
                brain_area = 'IPS';
            end
            [shuff]  = get_shuffled_activity(monkey,brain_area,p);

%%  plot divergence of seqPC1
            [~,reps] = size(divergepsth_seqPC1);
            ys = mean_divergepsth_seqPC1;
            % resign the sign to prefer direction
            abs_value = abs(ys);
            maxindex = find(abs_value == max(abs_value));
            if ys(maxindex)>0
                ys = ys;
                divergepsth_seqPC1 = divergepsth_seqPC1;
            else
                ys = -ys;
%                 divergepsth_seqPC1 = -divergepsth_seqPC1;
                divergepsth_seqPC1 = cellfun(@(x) -x, divergepsth_seqPC1,'UniformOutput', false);
            end 
            num_select = 200;
%             num_select = 15;
            select_divergepsth_seqPC1 = divergepsth_seqPC1(:, randperm(reps, num_select));
            divergepsth_seqPC1_temp =cellfun(@(x) mean(x,2), select_divergepsth_seqPC1,'UniformOutput', false);
            mat_divergepsth_seqPC1  = cell2mat(divergepsth_seqPC1_temp);  
            sem = nanstd(mat_divergepsth_seqPC1,0,2)./sqrt(num_select);
            CI = sem*1.96;
            
            % for shuffle part
%             [~,pp1] = cellfun(@(x) ttest2(x',shuff,'Tail','right'), select_divergepsth_seqPC1, 'UniformOutput', false);
%             [~,pp2] = cellfun(@(x) ttest2(x',-shuff,'Tail','left'), select_divergepsth_seqPC1, 'UniformOutput', false);
            
            base_line = mean(mean(shuff));
%             [~,pp1] = cellfun(@(x) ttest(x'-base_line,0,'Tail','right'), select_divergepsth_seqPC1, 'UniformOutput', false);
%             [~,pp2] = cellfun(@(x) ttest(x'-base_line,0,'Tail','left'), select_divergepsth_seqPC1, 'UniformOutput', false);
            [~,pp1] = cellfun(@(x) ttest(x',base_line,'Tail','right'), select_divergepsth_seqPC1, 'UniformOutput', false);
%             [~,pp2] = cellfun(@(x) ttest(x',-base_line,'Tail','left'), select_divergepsth_seqPC1, 'UniformOutput', false);
    
            
            [~,ppp1] = ttest2(mat_divergepsth_seqPC1', shuff,'Tail','right');
            [~,ppp2] = ttest2(mat_divergepsth_seqPC1', -shuff,'Tail','left');
            if ibrain_area==2
%                 temp = divergepsth_seqPC1>mean_shu_seqPC1;
%                 temp = cellfun(@(x,y) x<p_critical|y<p_critical, pp1, pp2,'UniformOutput', false);
                temp = cellfun(@(x,y) x<p_critical, pp1, 'UniformOutput', false);
%                 col_indices = cellfun(@(x) find(x, 1, 'first'), temp, 'UniformOutput', false);
                col_indices = cellfun(@(x) ...
                        find(conv(double(x), [1 1 1 1 1], 'valid') == 5, 1, 'first'), ...
                        temp, ...
                        'UniformOutput', false);

                MST_starttime = arrayfun(@(i) time(col_indices{i}),  1:length(col_indices), 'UniformOutput', false);
                
                fprintf('MST,  median=%d\n', median(cell2mat(MST_starttime)));


                filename = ['./timedata/',monkey,'p',num2str(p),'MSTtime.mat'];
                save(filename,"MST_starttime")
            elseif ibrain_area==4
%                 temp = divergepsth_seqPC1>mean_shu_seqPC1;
%                 temp = cellfun(@(x,y) x<p_critical|y<p_critical, pp1, pp2,'UniformOutput', false);
                temp = cellfun(@(x,y) x<p_critical, pp1, 'UniformOutput', false);
                col_indices = cellfun(@(x) find(x, 1, 'first'), temp, 'UniformOutput', false);
                LIP_starttime = arrayfun(@(i) time(col_indices{i}),  1:length(col_indices), 'UniformOutput', false);
                
                fprintf('LIP,  median=%d\n', median(cell2mat(LIP_starttime)));
         
                filename = ['./timedata/',monkey,'p',num2str(p),'LIPtime.mat'];
                save(filename,"LIP_starttime")
            end            
            significant_loc = (ppp1<p_critical)|(ppp2<p_critical);
            if ~isempty( LIP_starttime)
                if ~all(cellfun(@isempty, MST_starttime))&&~all(cellfun(@isempty, LIP_starttime))
%                     start_p = ranksum(cell2mat(MST_starttime),cell2mat(LIP_starttime));
                    [~,start_p] = ttest2(cell2mat(MST_starttime),cell2mat(LIP_starttime));

                else
                    start_p = nan;
                end
            end
            this_color = min([color(1,:)+[1,1,1]*1*(1-1*ibrain_area/4);1,1,1]);
            if ibrain_area==2
                if p==4
                    this_color = [30,114,255]/255;
                    MST_color = this_color;
                else
                    MST_color = this_color;
                end
                
            elseif ibrain_area==4
                if p==4
                    this_color = [30,64,139]/255;
                    LIP_color = this_color;
                else
                    LIP_color = this_color;
                end
            end

            % 箱线图（'Colors' 设置颜色）
            if ibrain_area==4                
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
                title_info = ['p=',num2str(start_p)];
                title(title_info)
                SetFigure
%                 str = sprintf('MST median = %d\n' , median(cell2mat(MST_starttime)));
%                 annotation('textbox', [0.15 0.85 0.35 0.058], 'String', str);
%                 str = sprintf('LIP median = %d\n' , median(cell2mat(LIP_starttime)));
%                 annotation('textbox', [0.15 0.55 0.35 0.058], 'String', str);
                str = sprintf('MST mean = %d\n' , round(mean(cell2mat(MST_starttime))));
                annotation('textbox', [0.15 0.55 0.35 0.058], 'String', str);
                str = sprintf('LIP mean = %d\n' , round(mean(cell2mat(LIP_starttime))));
                annotation('textbox', [0.15 0.85 0.35 0.058], 'String', str);
            end
            %{
            % bar 图
            if ibrain_area==4                
                figure(p*10+3)
                MST_starttime_val = cell2mat(MST_starttime);
                LIP_starttime_val = cell2mat(LIP_starttime);
               
               % 手动指定bin边缘（根据你的数据范围调整）
                all_data = [MST_starttime_val, LIP_starttime_val];
               data_range = max(all_data) - min(all_data);
               start_edge = floor(min(all_data));
               end_edge = ceil(max(all_data));
               
               % 创建固定宽度的bin边缘
               bin_width = 20; % 设置你想要的bin宽度
               edges = start_edge:bin_width:end_edge;
               
               hold on
               histogram(MST_starttime_val, edges, 'FaceColor', MST_color, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
               histogram(LIP_starttime_val, edges, 'FaceColor', LIP_color, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
               
               % 美化图形
               xlabel('Response');
               ylabel('Count');
               box off;
               set(gca, 'TickDir', 'out');
               title_info = ['p=', num2str(start_p)];
               title(title_info)
               legend({'Group 1', 'Group 2'});
               SetFigure
            end
            %}
            figure(p*10+1)
            hold on
            %     plot(time,ys,'Color',this_color,'LineWidth',linesize)
            if ~isnan(ys)
                shadedErrorBar(time,ys,CI,'lineprops',{'Color',this_color,'LineStyle',lineStyle,'LineWidth',linesize},'transparent',transparent);
                hold on
                shuffle_ys = nanmean(shuff);
                shuffle_CI = nanstd(shuff)./size(shuff,1)*1.96;
                shadedErrorBar(time,shuffle_ys,shuffle_CI,'lineprops',{'Color',this_color,'LineStyle','--'},'transparent',transparent);
%                 plot(time,ones(size(time))*mean(shuff),'Color',this_color)
%                 plot(time,-ones(size(time))*mean(shuff),'Color',this_color)
%                 plot(shu_time,shuffle_ys,'Color',this_color)
%                 plot(shu_time,-shuffle_ys,'Color',this_color)
%                 plot(time(significant_loc),zeros(size(time(significant_loc)))-ibrain_area*0.01,'*','Color',this_color)
                width =20; height = 0.004;
                plotsignificant(time(significant_loc),zeros(size(time(significant_loc)))-ibrain_area*0.01,width,height,this_color)
                title('seqPCAdim 1' )
            end
            ylim([-0.1,0.4])
            SetFigure            
        end
    end
end