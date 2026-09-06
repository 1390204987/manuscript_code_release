% plot single variable signed variable divergence
% by zn 2025/4/12
%% Load data
clear; close all;
%
monkey = 'F';
brain_area = 'IPS';

STS_data = load([monkey,'_',brain_area,'color_all_20ms']);
STS_seqPCA_vec = load([monkey,'_colorvariable_all_',brain_area,'.mat']);
STS_mTDR_vec = load([monkey,'_',brain_area,'color_all_mTDR.mat']);
% shu_STS_data = load('F_IPScolor_shuffleall_all.mat');
% shu_STS_seqPCA_vec = load('F_colorvariable_shuffleall_IPS.mat');
% shu_STS_mTDR_vec = load('F_IPScolor_shuffleall_mTDR.mat');
%
% IPS_data = load('D_IPScolor_all');
% IPS_seqPCA_vec = load('D_colorvariable_all_IPS.mat');
% IPS_mTDR_vec = load('D_IPScolor_all_mTDR.mat');
% shu_IPS_data = load('D_IPScolor_shuffleall_all.mat');
% shu_IPS_seqPCA_vec = load('D_colorvariable_shuffleall_IPS.mat');
% shu_IPS_mTDR_vec = load('D_IPScolor_shuffleall_mTDR.mat');
color_{1,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
color_{2,1} = [0,191,255;155,48,255]/255; %color
color_{3,1} = [139,10,80;72,61,139]/255; %choice
color_{4,1} = [0,0,139;205,133,0]/255; %saccade
color_{5,1} = [0,0,139;205,133,0]/255; %last saccade
color_{6,1} = [0,0,1;0,0,0.9;0,0,0.8;0,0,0.7;0,0,0.6]; %coherence
color_{7,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
linesize = 2;
markersize = 10;
transparent = 1;
lineStyle = '-';
variable_name={'heading';'color';'choice';'saccade';'last_sac';'coherence'};
p_critical = 0.001;
for p = 6%1:4
    
    color = color_{p,1};

    data =STS_data; vec_seqPCA = STS_seqPCA_vec; vec_mTDR = STS_mTDR_vec;
%     shu_data =shu_STS_data; shu_vec_seqPCA = shu_STS_seqPCA_vec; shu_vec_mTDR = shu_STS_mTDR_vec;
%     data =IPS_data; vec_seqPCA = IPS_seqPCA_vec; vec_mTDR = IPS_mTDR_vec;
%     shu_data =shu_IPS_data; shu_vec_seqPCA = shu_IPS_seqPCA_vec; shu_vec_mTDR = shu_IPS_mTDR_vec;
    
    not_nan_columns = all(~isnan(vec_seqPCA.variablePC{1,p}), 1);
    vec_dims = sum(not_nan_columns);
    real_dim{p} = []; % if a dim also has significant difference with shuffle, it will be set to real dim
    mean_divergencepsth_seqPC = [];
    if ~isempty(data)        
        time = data.t_centers{:};
        if p == 3
            [psth,heading_cons,choice_cons,unique_cohe] = get_2v_PSTH(data.Xmat,data.Ymat,data.hk);
%             [shu_psth,heading_cons,choice_cons,unique_cohe] = get_2v_PSTH(shu_data.Xmat,shu_data.Ymat,shu_data.hk);
        else
            [psth,~,~,unique_cons,unique_cohe] = get_PSTH(data.Xmat,data.Ymat,p,data.hk);
%             [shu_psth,~,~,unique_cons,unique_cohe] = get_PSTH(shu_data.Xmat,shu_data.Ymat,p,shu_data.hk);
        end

        coherence_0=0;
        for idim = 1:vec_dims
            if p == 3
                [mean_divergepsth_seqPC{idim,1},divergepsth_seqPC{idim,1}]=getdivergence_new(time,psth,choice_cons,heading_cons,unique_cohe,vec_seqPCA.variablePC{1,p}(:,idim),coherence_0,color,linesize);
%                 [mean_shu_divergepsth_seqPC1,shu_divergepsth_seqPC1]=getdivergence_new(time,shu_psth,choice_cons,heading_cons,unique_cohe,shu_vec_seqPCA.variablePC{1,p}(:,1),coherence_0,color,linesize);
                
            else
                [mean_divergepsth_seqPC{idim,1},divergepsth_seqPC{idim,1}]=getdivergence_new(time,psth,unique_cons,[],unique_cohe,vec_seqPCA.variablePC{1,p}(:,idim),coherence_0,color,linesize);
%                 [mean_shu_divergepsth_seqPC1,shu_divergepsth_seqPC1]=getdivergence_new(time,shu_psth,unique_cons,[],unique_cohe,shu_vec_seqPCA.variablePC{1,p}(:,1),coherence_0,color,linesize);
            end
            [shuff]  = get_shuffled_activity(monkey,brain_area,p);
            [~,reps] = size(divergepsth_seqPC{idim,1});
            ys = mean_divergepsth_seqPC{idim};   
            
            stimulus_period = time>0&time<1500;
            ys_stimulus_period = ys(stimulus_period);
            abs_value = abs(ys_stimulus_period);
            maxindex = find(abs_value == max(abs_value));
            if ys_stimulus_period(maxindex)>0
                ys = ys;
                divergepsth_seqPC{idim} = divergepsth_seqPC{idim};
                variablePC{1,p}(:,idim) = vec_seqPCA.variablePC{1,p}(:,idim);
            else
                ys = ys;
%                 divergepsth_seqPC{idim} = -divergepsth_seqPC{idim};
                divergepsth_seqPC{idim} = cellfun(@(x) -x, divergepsth_seqPC{idim},'UniformOutput',false);
                variablePC{1,p}(:,idim) = -vec_seqPCA.variablePC{1,p}(:,idim);
            end
            temp = cellfun(@(x) mean(x,2),divergepsth_seqPC{idim},'UniformOutput', false);
            mat_divergepsth_seqPC{idim} = cell2mat(temp);
            sem = nanstd( mat_divergepsth_seqPC{idim},0,2)./size(mat_divergepsth_seqPC{idim},2);
%             sem = nanstd(divergepsth_seqPC{idim},0,2)./size(mat_divergepsth_seqPC{idim},2);
            CI = sem;
%             shuffle_ys = mean_shu_divergepsth_seqPC1;
%             shu_abs_value = abs(shuffle_ys);
%             shu_maxindex = find(shu_abs_value == max(shu_abs_value));
%             if shuffle_ys(shu_maxindex)>0
%                 shuffle_ys = shuffle_ys;
%                 shu_divergepsth_seqPC1 = shu_divergepsth_seqPC1;
%             else
%                 shuffle_ys=-shuffle_ys;
%                 shu_divergepsth_seqPC1 = -shu_divergepsth_seqPC1;
%             end
            % 使用 arrayfun 生成 1000 次随机打乱
%             num_shuffles = 100;
%             shuffled_indices = arrayfun(@(x) randperm(length(time)), 1:num_shuffles, 'UniformOutput', false);
%             temp_shuffled = cellfun(@(idx) shuffle_ys(idx)', shuffled_indices, 'UniformOutput', false);
%             shu_shu_seqPC1 = cell2mat(temp_shuffled);
            base_line = mean(mean(shuff));
            [~,pp1] = ttest(mat_divergepsth_seqPC{idim}',base_line,'Tail','right');
            [~,pp2] = ttest(mat_divergepsth_seqPC{idim}',-base_line,'Tail','left');
            
            pp = (pp1<p_critical)|(pp2<p_critical);
            significant_loc = (pp1<p_critical)|(pp2<p_critical);
            this_color = min([color(1,:)+[1,1,1]*0.75*(1-1*1);1,1,1]);
            
            if ~isnan(ys)
                figure()
                shadedErrorBar(time,ys,CI,'lineprops',{'Color',this_color,'LineStyle',lineStyle},'transparent',transparent);
                hold on
%                 shadedErrorBar(time,shuffle_ys,shuffle_CI,'lineprops',{'Color',this_color,'LineStyle','--'},'transparent',transparent);
%                 shadedErrorBar(time,-shuffle_ys,shuffle_CI,'lineprops',{'Color',this_color,'LineStyle','--'},'transparent',transparent);
                plot(time,mean(shuff,1))
                plot(time,-mean(shuff,1))
                plot(time(significant_loc),zeros(size(time(significant_loc))),'*','Color',this_color)
                SetFigure
                titlename = ['seqPCA in dim',num2str(idim)];
                title(titlename)
            end
            
            if sum(pp(stimulus_period))>4
                real_dim{p} = [real_dim{p},idim];
            end
            
        end


    end
end
filename = 'F_colorvariable_all_IPS_redit.mat';
variablePC = variablePC;cellID_pick = vec_seqPCA.cellID_pick;
save(filename,'variablePC','cellID_pick','real_dim')
