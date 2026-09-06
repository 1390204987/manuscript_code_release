% plot two monkey data on same figure
%%
clear; close all;
for M_data = 1:2
    if M_data==1
%         FolderPath_HD={'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorIPS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorIPS_NP\'};
        FolderPath_HD ={'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorSTS_NP\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorNTSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTSTS_NP\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorNTIPS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTIPS_NP\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorMST_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorMST_NP\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorNTMST_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTMST_NP\'};
        FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dIPS_LA\'};
    else
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseIPS_LA\';};
        FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseSTS_LA\';};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseVIP_LA\';};
%         FolderPath_HD ={'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorSTS_NP\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorNTSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTSTS_NP\'};
        FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dIPS_LA\'};
    end
    
    if length(FolderPath_HD)==1
        Files_HD = dir(fullfile(FolderPath_HD{1},'*PSTH.mat'));
        pathfilenum = length(Files_HD);
    else
        Files_HD = []; pathfilenum=[];
        for ipath = 1:length(FolderPath_HD)
            ipathfiles = dir(fullfile(FolderPath_HD{ipath},'*PSTH.mat'));
            Files_HD = [Files_HD;ipathfiles];
            pathfilenum(ipath) = length(ipathfiles);
        end
        pathfilenum = cumsum(pathfilenum);
    end
    File_work_count = 0;
     list=[];cellnum=[];blocknum=[];cellID=[];
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
        monkeyID = str2num(list{1,ct}.result.FILE(charloc(1)+1:charloc(2)-1));
        cellnum(ct) = str2double(list{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;
        blocknum(ct) = str2double(list{1,ct}.result.FILE(charloc(3)+1:end));
        cellID(ct) = str2double(list{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;   %+blocknum(ct)/10;
    end
    %%
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
     list_1D=[];p_1D=[];az_1D=[];cellnum_1D=[];cellID_1D=[];
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
    list_1D = list_1D(p_1D<1);
    
 %%   
%     [cellID_pick,cellindex,cellindex_1D] = intersect(cellID,select_cellID_1D);
    [cellID_pick,cellindex,cellindex] = intersect(cellID,cellID);
    
    %% pack data
     heading_per_trial=[]; coherence_per_trial=[]; relativeheading_per_trial=[]; choice_per_trial=[]; choice_last_trial=[]; targetloc_per_trial=[];
      signtargetloc_per_trial=[]; signtargetloc_last_trial=[]; saccade_per_trial=[]; spike_hist=[]; spike_aligned=[]; spike_per_cell=[];
    for f = 1:length(cellID_pick)
        heading_per_trial{f,1} = (list{1,cellindex(f)}.result.heading_per_trial)';
        if monkeyID~=8 % 8 is yuri
            coherence_per_trial{f,1} = (list{1, cellindex(f)}.result.coherence_per_trial)';
        else
            coherence_per_trial(f,1) = ones(size(heading_per_trial{f,1}));
        end
        if length(unique(coherence_per_trial{f,1})) >1&monkeyID~=8
            relativeheading_per_trial{f,1} = sign(heading_per_trial{f,1}).*coherence_per_trial{f,1};
        else
            relativeheading_per_trial{f,1} = heading_per_trial{f,1};
        end
        choice_per_trial{f,1} = (list{1, cellindex(f)}.result.choice_per_trial)';
        choice_last_trial{f,1} = [0;0;choice_per_trial{f,1}(2:end-1)];
        targetloc_per_trial{f,1} = (list{1, cellindex(f)}.result.targetloc_per_trial)';
        signtargetloc_per_trial{f,1} = sign(targetloc_per_trial{f,1});
        signtargetloc_last_trial{f,1} = [0;0;signtargetloc_per_trial{f,1}(2:end-1)];
    %     signtargetloc_per_trial{f,1} = Shuffle(signtargetloc_per_trial{f,1});
        saccade_per_trial{f,1} = sign(list{1, cellindex(f)}.result.targetloc_per_trial.*list{1, cellindex(f)}.result.choice_per_trial)';

        spike_hist{f,1} = [{list{1,cellindex(f)}.result.spike_hist{1,1}},list{1,cellindex(f)}.result.t_centers{1,1}'];
        spike_aligned{f,1} = [{list{1,cellindex(f)}.result.spike_aligned{1,1}},{list{1,cellindex(f)}.result.spike_aligned{2,1}}];

        %for time window selection
        stim_dur(f,1) = list{1,cellindex(f)}.result.align_offsets_others{1,2}(1,2) - list{1,cellindex(f)}.result.align_offsets_others{1,2}(1,1);
        stimulus_on(f,1) = 0;
        stimulus_end(f,1) = 0+stim_dur(f,1);

        %method1 calcalate with psth data
        time_window1(f,1) = find(spike_hist{f,1}{1,2}(:) > stimulus_on(f,1)+0,1);
        time_window1(f,2) = find(spike_hist{f,1}{1,2}(:) < stimulus_end(f,1),1,'last');
        %method2 calculate with 0/1 spike data
        time_window2(f,1) = find(spike_aligned{f,1}{1,2}(:) > stimulus_on(f,1)+500,1);
        time_window2(f,2) = find(spike_aligned{f,1}{1,2}(:) < stimulus_end(f,1),1,'last');

        spike_per_cell{f,1} = mean(spike_hist{f,1}{1,1}(:,time_window1(f,1):time_window1(f,2)),2);
    end
    
    sensory_period = [0,1500];
    color_period = [0,1500];
    spacial_period = [0,1500];
    choice_period = [0,1500];
    
    time_window_forheading(1,1) = find(spike_hist{1,1}{1,2}(:) < stimulus_on(1,1)+sensory_period(1),1);
    time_window_forheading(1,2) = find(spike_hist{1,1}{1,2}(:) < stimulus_on(1,1)+sensory_period(2),1,'last');
    
    time_window_forcolor(1,1) = find(spike_hist{1,1}{1,2}(:) < stimulus_on(1,1)+color_period(1),1);
    time_window_forcolor(1,2) = find(spike_hist{1,1}{1,2}(:) < stimulus_on(1,1)+color_period(2),1,'last');
    
    time_window_forspacial(1,1) = find(spike_hist{1,1}{1,2}(:) < stimulus_on(1,1)+spacial_period(1),1);
    time_window_forspacial(1,2) = find(spike_hist{1,1}{1,2}(:) < stimulus_on(1,1)+spacial_period(2),1,'last');
    
    time_window_forchoice(1,1) = find(spike_hist{1,1}{1,2}(:) < stimulus_on(1,1)+choice_period(1),1);
    time_window_forchoice(1,2) = find(spike_hist{1,1}{1,2}(:) < stimulus_on(1,1)+choice_period(2),1,'last');

    %% 2.0 choice probability
    LEFT = -1;
    RIGHT = 1;
    spike_per_cell_forheading=[]; spike_per_cell_forcolor=[]; spike_per_cell_forspacial=[]; spike_per_cell_forchoice=[];
    zleftheading_spike_per_cell=[]; zrightheading_spike_per_cell=[]; roc_heading=[];  line=[]; T1select=[];
    resp_left_all=[]; resp_right_all=[]; Zchoicespike_per_cell=[];  resp_left_choose=[];  resp_right_choose=[]; CP=[];  CP_all=[];
    choice_roc_T1L=[]; choice_roc_T1R=[]; CP_T1L=[];  CP_T1R=[]; CP_T1L_true=[]; CP_T1R_true=[]; abs_CP_diverge=[]; pvalue_CP_diverge=[];
    z_firing_rates=[]; firing_rates=[];sacROC=[];
    for f = 1:length(cellID_pick)
%         spike_per_cell_forheading{f,1} = mean(spike_hist{f,1}{1,1}(:,time_window_forheading(1,1):time_window_forheading(1,2)),2);
%         spike_per_cell_forcolor{f,1} = mean(spike_hist{f,1}{1,1}(:,time_window_forcolor(1,1):time_window_forcolor(1,2)),2);
%         spike_per_cell_forspacial{f,1} = mean(spike_hist{f,1}{1,1}(:,time_window_forspacial(1,1):time_window_forspacial(1,2)),2);
%         spike_per_cell_forchoice{f,1} = mean(spike_hist{f,1}{1,1}(:,time_window_forchoice(1,1):time_window_forchoice(1,2)),2);
        
        firing_rates{f,1} = sum(spike_aligned{f,1}{1,1}(:,time_window2(1,1):time_window2(1,2)),2)/((time_window2(1,2)-time_window2(1,1))/1000);

        [rr,pp] = corrcoef(relativeheading_per_trial{f,1}(:), firing_rates{f,1}(:));
        
        line(f,1) = rr(1,2);
        line(f,2) = pp(1,2);
        % roc for heading-------------------------------------------------------
        
        % roc for choice-----------------------------------------------------------
        see_targetloc = signtargetloc_per_trial{f,1};
        %    see_targetloc = signtargetloc_last_trial{f,1};
        for loc = 1:3
            if loc == 1
                T1select{f,loc} = (see_targetloc == -1); %red on left
            elseif loc == 2
                T1select{f,loc} = (see_targetloc == 1);
            else
                T1select{f,loc} = (see_targetloc == -1|see_targetloc == 1);
            end
            
            unique_relativeheading = unique(relativeheading_per_trial{f,1});
            see_choice = choice_per_trial{f,1};
            %        see_choice = choice_last_trial{f,1};
            %        see_choice = choice_last_trial{f,1};
            for i = 1:length(unique_relativeheading)
                select = (relativeheading_per_trial{f,1} == unique_relativeheading(i));
                z_choice = firing_rates{f,1}(select);
                z_choice = (z_choice - mean(z_choice))/std(z_choice);
                z_firing_rates{f,1}(select,1) = z_choice;
                
                resp_left_choose{loc}{f,i} = firing_rates{f,1}(select & (see_choice == LEFT) & T1select{f,loc} );
                resp_right_choose{loc}{f,i} = firing_rates{f,1}(select & (see_choice == RIGHT) & T1select{f,loc} );
                
                if (length(resp_left_choose{loc}{f,i}) <= 3) || (length(resp_right_choose{loc}{f,i}) <= 3)   % make sure each stim_type has at least 3 data values
                    resp_left_choose{loc}{f,i}(select) = 9999;   % similar to NaN, just make a mark
                    resp_right_choose{loc}{f,i}(select) = 9999;
                end
                
                % Each unique heading
                if (length(resp_left_choose{loc}{f,i}) > 3) && (length(resp_right_choose{loc}{f,i}) > 3)
                    [CP{loc,i}(f,1),~ ,perm] = rocN( resp_left_choose{loc}{f,i},resp_right_choose{loc}{f,i});  % raw firing rate
                    [~,CP{loc,i}(f,2)] = ttest2(resp_left_choose{loc}{f,i},resp_right_choose{loc}{f,i});
                    %                 CP{loc,i}(f,2) = perm.pValue;
                else
                    CP{loc,i}(f,1) = NaN;
                    CP{loc,i}(f,2) = NaN;
                end
                
                if logical(line(f,1)>0)&&logical(~isnan(CP{loc,i}(f,1)))
                    %             if logical(roc_heading(f,1)<0.5)&&logical(~isnan(CP{loc,i}(f,1)))
                    CP{loc,i}(f,1) = 1-CP{loc,i}(f,1);
                end
                
            end
            % now across all data (z-score for grand CP? HH)
            % Grand CP (grand cp not conditioned on specific heading thus the value might be affected by heading selectivity ?)
            resp_left_all{f,loc} = z_firing_rates{f,1}((see_choice == LEFT) & (z_firing_rates{f,1}~=9999)& T1select{f,loc} );
            resp_right_all{f,loc} = z_firing_rates{f,1}((see_choice == RIGHT) & (z_firing_rates{f,1}~=9999)& T1select{f,loc}  );
            
            
            if (length(resp_left_all{f,loc}) > 3) && (length(resp_right_all{f,loc}) > 3)
                [CP_all{loc,1}(f,1),~ ,perm] = rocN( resp_left_all{f,loc},resp_right_all{f,loc},100,1000);   % z-scored firing rate
                CP_all{loc,1}(f,2) = perm.pValue;
%                 [~,CP_all{loc,1}(f,2)] = ttest2(resp_left_all{f,loc},resp_right_all{f,loc});
            else
                CP_all{loc,1}(f,1) = NaN;
                CP_all{loc,1}(f,2) = NaN;
            end
            
            roc_choice{loc,1}(f,1) = CP_all{loc,1}(f,1);
            roc_choice{loc,1}(f,2) = CP_all{loc,1}(f,2);
            
            if  logical(line(f,1)>0)&&logical(~isnan(CP_all{loc,1}(f,1)))
                %        if  logical(roc_heading(f,1)<0.5)&&logical(~isnan(CP_all{loc,1}(f,1)))
                CP_all{loc,1}(f,1) = 1-CP_all{loc,1}(f,1);
            end
        end
        %{
        num_shuffles = 1000;
        num_trials = length(choice_per_trial{f,1}); % 获取 choice_per_trial 的长度
        % 使用 arrayfun 生成 1000 次随机打乱
        shuffled_indices = arrayfun(@(x) randperm(num_trials), 1:num_shuffles, 'UniformOutput', false);
        % 如果需要打乱后的 choice_per_trial 数据
        shuffled_T1loc = cellfun(@(idx)  signtargetloc_per_trial{f,1}(idx), shuffled_indices, 'UniformOutput', false);
        T1loc = { signtargetloc_per_trial{f,1},shuffled_T1loc{:}}; % the first column is true lable the others is shuffled label
        Zresp_left_choice_T1L = cellfun(@(x) Zchoicespike_per_cell{f,1}((choice_per_trial{f,1}==-1)&(sign(x)==-1),:),T1loc, 'UniformOutput', false);
        Zresp_right_choice_T1L = cellfun(@(x) Zchoicespike_per_cell{f,1}((choice_per_trial{f,1}==1)&(sign(x)==-1),:),T1loc,'UniformOutput', false);
        
        Zresp_left_choice_T1R = cellfun(@(x) Zchoicespike_per_cell{f,1}((choice_per_trial{f,1}==-1)&(sign(x)==1),:),T1loc,'UniformOutput', false);
        Zresp_right_choice_T1R = cellfun(@(x) Zchoicespike_per_cell{f,1}((choice_per_trial{f,1}==1)&(sign(x)==1),:),T1loc,'UniformOutput', false);
        
        
        choice_roc_T1L{f,1}= cellfun(@(x,y) rocN(x,y),Zresp_left_choice_T1L,Zresp_right_choice_T1L);  %set the left as x axis
        choice_roc_T1R{f,1}= cellfun(@(x,y) rocN(x,y),Zresp_left_choice_T1R,Zresp_right_choice_T1R);  %set the left as x axis
        
        if  logical(line(f,1)>0)%&logical(~isnan(choice_roc_T1L))
            CP_T1L{f,1} = 1-choice_roc_T1L{f,1};
            CP_T1R{f,1} = 1-choice_roc_T1R{f,1};
        else
            CP_T1L{f,1} = choice_roc_T1L{f,1};
            CP_T1R{f,1} = choice_roc_T1R{f,1};
        end
        CP_T1L_true(f,1) = CP_T1L{f,1}(1,1);
        CP_T1R_true(f,1) = CP_T1R{f,1}(1,1);
        
        abs_CP_diverge{f,1} = abs(CP_T1L{f,1}-CP_T1R{f,1});
        pvalue_CP_diverge(f,1) = sum(abs_CP_diverge{f,1}(:,2:end)>abs_CP_diverge{f,1}(:,1),2)/num_shuffles;
        %}
        
        %calculate saccade ROC
        see_sac = saccade_per_trial{f,1};
        resp_leftsac_all{f,1} = z_firing_rates{f,1}((see_sac == LEFT) & (z_firing_rates{f,1}~=9999));
        resp_rightsac_all{f,1} = z_firing_rates{f,1}((see_sac == RIGHT) & (z_firing_rates{f,1}~=9999));
        
        if (length(resp_leftsac_all{f,1}) > 3) && (length(resp_rightsac_all{f,1}) > 3)
            [sacROC(f,1),~ ,perm] = rocN( resp_leftsac_all{f,1},resp_rightsac_all{f,1},100,1000);   % z-scored firing rate
            sacROC(f,2) = perm.pValue;
%             [~,sacROC(f,2)] = ttest2(resp_leftsac_all{f,1},resp_rightsac_all{f,1});
        else
            sacROC(f,1) = NaN;
            sacROC(f,2) = NaN;
        end
        
    end
    %%  3.0 plot figure
    % 3.2 plot choice probability
    showglobalCP = 1;
    if showglobalCP ==1
        for igroup = 1:length(CP_all)
            %         CP_plot{igroup,1} = CP_all{igroup,1}(line(:,2)<0.05,:);
            CP_plot{igroup,1} = CP_all{igroup,1}(:,:);
        end
    else
        CP_plot = CP(:,5);
    end
    %%
    % 3.3 plot cell's CP correlation, which calculated on different target location
    if length(unique(targetloc_per_trial{1,1}))>1
        if M_data==1 %monkey Drogon
            color = [30 144 255]/255;
            text_yloc = 0.8;
        elseif M_data==2 %monkey Friday
            color = [255 165 0]/255;
            text_yloc = 0.2;
        end
        figure(103)
        hold on
        if sum((CP_plot{1,1}(:,2) > 0.05) & (CP_plot{2,1}(:,2) > 0.05)) ~= 0
            plot_roc_cor(CP_plot{1,1}(:,1),CP_plot{2,1}(:,1), 1,text_yloc,{'o'},{color});
%             plot(CP_plot{1,1}(pvalue_CP_diverge<0.05,1),CP_plot{2,1}(pvalue_CP_diverge<0.05,1), 'o','MarkerEdgeColor',color,'MarkerFaceColor',color);        
            plot(CP_plot{1,1}(sacROC(:,2)<0.05,1),CP_plot{2,1}(sacROC(:,2)<0.05,1), 'o','MarkerEdgeColor',color,'MarkerFaceColor',color);       
        end
        SetFigure
    end
    fprintf('Monkey=%d,  roc_num=%d\n', M_data , sum(sacROC(:,2)<0.05));
%     significant_num = sum(pvalue_CP_diverge<0.05)
%     significant_percent = sum(pvalue_CP_diverge<0.05)/length(pvalue_CP_diverge)
    %% marker specific cell in the figure
%     cell1 = 11000018; %m24c115r3_293
%     cell2 = 11000023;%m24c110r3_383
    cell1 = 11500238; %m24c115r3_293
    cell2 = 11500293;%m24c110r3_383
    cell1_index = find(cellID_pick==cell1);
    cell2_index = find(cellID_pick==cell2);
    if ~isempty(cell1_index)
        cell1_x = CP_plot{1,1}(cell1_index,1);
        cell1_y = CP_plot{2,1}(cell1_index,1);
    end
    if ~isempty(cell2_index)
        cell2_x = CP_plot{1,1}(cell2_index,1);
        cell2_y = CP_plot{2,1}(cell2_index,1);
    end
    figure(103); hold on;
    pbaspect([1 1 1])
    xlim([0,1])
    yticks([0.1,0.5,0.9])
    xticks([0.1,0.5,0.9])
%     plot(cell1_x,cell1_y, 'v','MarkerEdgeColor',[0 0 0]);
%     plot(cell2_x,cell2_y, 's','MarkerSize',8,'MarkerEdgeColor',[0 0 0]);
%     plot(cell1_x,cell1_y, 'd','MarkerEdgeColor',[0 0 0]);
%     plot(cell2_x,cell2_y, 'h','MarkerSize',8,'MarkerEdgeColor',[0 0 0]);
end

