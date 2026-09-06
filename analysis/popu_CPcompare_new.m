,%2024/9/10 population CP Dynamic 
% function popu_CP(varargin)
% if ~isempty(varargin)
%     list = varargin{1};
% else
 close all;
 clear;
for iarea = 1:2
%     FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorIPS_LA\';'Z:\Data\Tempo\Batch\fingerpost\DrogoncolorIPS\'};
%     FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorSTS_LA\';'Z:\Data\Tempo\Batch\fingerpost\DrogoncolorMST\';'Z:\Data\Tempo\Batch\fingerpost\DrogoncoarsecolorMST\'};
%     FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorSTS_LA\'};
%     FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorIPS_LA\'};
    FolderPath_1D = [];
    FolderPath_HD= [];
    if iarea == 1
        FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\Friday1dSTS_NP\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorSTS_NP\'};
        FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorNTSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTSTS_NP\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_NP\simultaneous\MST_TOE\'};
%         FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dSTS_LA\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseSTS_LA\'};
    else
        FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dIPS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\Friday1dIPS_NP\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorIPS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorIPS_NP\'};
        FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorNTIPS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTIPS_NP\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_NP\simultaneous\LIP_TOE\'};
%         FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dIPS_LA\'};
%         FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseIPS_LA\'};
    end
    %% select cell significant in 1D tunning
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
    list = []; cellnum=[]; cellID=[];
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
        cellnum(ct) = str2double(list{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1));
        cellID(ct) = str2double(list{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;
    end
    
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
    list_1D=[]; p_1D=[]; az_1D=[]; cellnum_1D=[]; cellID_1D=[];
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
    select_cellID_1D = cellID_1D(p_1D<0.05);
%     [cellID_pick,cellindex,cellindex_1D] = intersect(cellID,select_cellID_1D);
    [cellID_pick,cellindex,cellindex] = intersect(cellID,cellID);
    length(cellID_pick)
%%
    heading_per_trial=[];coherence_per_trial=[];T1loc_per_trial=[];T1loc_last_trial=[];choice_per_trial=[];
    choice_last_trial=[];spike_hist=[]; spike_aligned=[];align_offsets_others=[];spike_rates=[];ts_per_trial=[];
    saccade_per_trial = [];
    for f = 1:length(cellID_pick)
        heading_per_trial{f,1} = (list{1,cellindex(f)}.result.heading_per_trial)';
        coherence_per_trial{f,1} = (list{1,cellindex(f)}.result.coherence_per_trial)';
        T1loc_per_trial{f,1} = (list{1,cellindex(f)}.result.targetloc_per_trial)';
        T1loc_last_trial{f,1} = [0;0;0;T1loc_per_trial{f,1}(1:end-3)];
        choice_per_trial{f,1} = (list{1,cellindex(f)}.result.choice_per_trial)';
        choice_last_trial{f,1} = [0;0;0;choice_per_trial{f,1}(1:end-3)];
        spike_hist{f,1} = list{1,cellindex(f)}.result.spike_hist;
        spike_aligned{f,1} = list{1,cellindex(f)}.result.spike_aligned;
        align_offsets_others{f,1} = list{1,cellindex(f)}.result.align_offsets_others;
        spike_rates{f,1} = list{1,cellindex(f)}.result.spike_rates;
        ts_per_trial{f,1} = list{1,cellindex(f)}.result.t_centers;
        saccade_per_trial{f,1} = sign((list{1,cellindex(f)}.result.targetloc_per_trial.*list{1,cellindex(f)}.result.choice_per_trial)');
    end
    align_markers = list{1,1}.result.align_markers;
    temp_duration_ratio = zeros(1,size(align_markers,1));
    for j = 1:size(align_markers,1)
    %     temp_duration_ratio(j) = align_markers{j,3} - align_markers{j,2};
            temp_duration_ratio(1) = 1600-(-200);
            temp_duration_ratio(2) = 200-(-300);
    end
    for ialign = 1:size(align_markers,1)
        if ialign==1
            %time window stimulus on aligned
            stim_dur = align_offsets_others{1,1}{1,2}(1,2) - align_offsets_others{1,1}{1,2}(1,1);
            stimulus_on = 0;
            stimulus_end = 0+stim_dur;
            %method1 calcalate with psth data
            time_window1(1,1) = find(ts_per_trial{1,1}{1,1}(:) > stimulus_on-800,1);
            time_window1(1,2) = find(ts_per_trial{1,1}{1,1}(:) < stimulus_end,1,'last');
            time_window1(1,3) = size(ts_per_trial{1,1}{1,1},2);
            %method2 calculate with 0/1 spike data
            time_window2(1,1) = find(spike_aligned{1,1}{2,1}(:) > stimulus_on+500,1);
            time_window2(1,2) = find(spike_aligned{1,1}{2,1}(:) < stimulus_end,1,'last');
            time_window2(1,3) = size(spike_aligned{1,1}{2,1},2);
        else
            % saccade on aligned
            %method1 calcalate with psth data
            time_window1(1,1) = find(ts_per_trial{1,1}{1,2}(:) > 0-600,1);
            time_window1(1,2) = find(ts_per_trial{1,1}{1,2}(:) < 0+500,1,'last');
            time_window1(1,3) = size(ts_per_trial{1,1}{1,2},2);
            %method2 calculate with 0/1 spike data
            time_window2(1,1) = find(spike_aligned{1,1}{2,2}(:) > 0-500,1);
            time_window2(1,2) = find(spike_aligned{1,1}{2,2}(:) < 0+500,1,'last');
            time_window2(1,3) = size(spike_aligned{1,1}{2,2},2);
        end
        relativeheading_per_trial=[];spikehist_per_trial=[];
        CP_T1L=[];CP_T1R=[];CP=[];
        for f = 1:length(cellID_pick)
            z_fir_per_trial=[];
            if length(unique(coherence_per_trial{f,1}))==1
                relativeheading_per_trial{f,1} = heading_per_trial{f,1};
            else
                relativeheading_per_trial{f,1} = sign(heading_per_trial{f,1}).*coherence_per_trial{f,1};
            end
            unirelativeheading = unique(relativeheading_per_trial{f,1});
            unichoice = unique(choice_per_trial{f,1});
            spikehist_per_trial{f,1} = spike_hist{f,1}{1,ialign};
            
            % calculate firing rate for select neuron

            fir_per_trial = sum(spike_aligned{f,1}{1,1}(:,time_window2(1,1):time_window2(1,2)),2)/((time_window2(1,2)-time_window2(1,1))/1000);

            %%
            % zscore
            % shift window, currently using 100ms shift with 1000ms width
            window_interval = 10;% 10*10ms
            step = 10; %10*10ms
            endloop = floor((time_window1(1,3)-time_window1(1,1)-window_interval/2)/step);
            t_center_step{1,ialign} = ts_per_trial{1,1}{1,ialign}(time_window1(1,1):step:time_window1(1,1)+step*endloop);
            Zspike_heading_pertrial = [];
            for iheading = 1:length(unirelativeheading)
                select = (relativeheading_per_trial{f,1} == unirelativeheading(iheading));
                z_heading_dist = spikehist_per_trial{f,1}(select,:);
                z_fir_per_trial(select,:) = (fir_per_trial(select,:)-mean(fir_per_trial(select,:)))/std(fir_per_trial(select,:));

                for s = 1:endloop+1
                    vertical = time_window1(1,1)-0.5*window_interval+step*(s-1): time_window1(1,1)+ 0.5*window_interval+step*(s-1) ;
                    if mean(z_heading_dist(:,vertical),2)==0
                        trial_num = size(z_heading_dist(:,vertical),1);
                        z_heading_dist_step(:,s) =zeros(trial_num,1);
                    else
                        z_heading_dist_step(:,s) = (mean(z_heading_dist(:,vertical),2) - mean(mean(z_heading_dist(:,vertical),2)))/std(mean(z_heading_dist(:,vertical),2),1);
                    end
                end
                Zspike_heading_pertrial(select,:) = z_heading_dist_step;
                z_heading_dist_step = [];
            end
            
                    see_T1loc = T1loc_per_trial{f,1};
%             see_T1loc = T1loc_last_trial{f,1};
                    see_choice = choice_per_trial{f,1};
%             see_choice = choice_last_trial{f,1};
                    see_sac = saccade_per_trial{f,1};
%                     shuffle_seesac = Shuffle(saccade_per_trial{f,1});
%                     shuffle_seechoice = Shuffle(see_choice);
            
            Zresp_left_choice = Zspike_heading_pertrial((see_choice==-1),:);
            Zresp_right_choice = Zspike_heading_pertrial((see_choice==1),:);
            
            arranged_choice = [see_choice(see_choice==-1,:);see_choice((see_choice==1),:)];
            arranged_choice_index = [find(see_choice==-1);find(see_choice==1)];
            
            
            if ialign ==1
                allZresp_left_choice = z_fir_per_trial((see_choice==-1),:);
                allZresp_right_choice = z_fir_per_trial((see_choice==1),:);
            end
            
            
            Zresp_left_sac = Zspike_heading_pertrial((see_sac==-1),:);
            Zresp_right_sac = Zspike_heading_pertrial((see_sac==1),:);
            
            arranged_sac = [see_sac(see_sac==-1,:);see_sac((see_sac==1),:)];
            arranged_sac_index = [find(see_sac==-1);find(see_sac==1)];
            if ialign ==1
                allZresp_left_sac = z_fir_per_trial((see_sac==-1),:);
                allZresp_right_sac = z_fir_per_trial((see_sac==1),:);
            end
            
            %% T test evaluate choice selectivity
            choice_roc_T1L = [];
            choice_roc_T1R = [];
            choice_roc = [];
            sac_roc = [];
            shuffle_sac_roc = [];
            shuffle_choice_roc = [];
            
            for s = 1:endloop+1

                if length(Zresp_left_choice(:,s)) >3 && length(Zresp_right_choice(:,s))>3
                    [choice_roc(s,1),~,perm] = rocN(Zresp_left_choice(:,s),Zresp_right_choice(:,s),100,1000);  %set the left as x axis
%                     [~,choice_roc(s,2)] = ttest2(Zresp_left_choice(:,s),Zresp_right_choice(:,s)); % Use ttest2 instead of permutation for speed (the outcomes are quite similar). HH20140601
                    choice_roc(s,2) = perm.pValue;
                    if ~isnan(perm.pValue)
                        shuffle_choice_roc(s,:) = perm.auROCPerm;
                    else
                        shuffle_choice_roc(s,:) = NaN(1000,1);
                    end
                else
                    [choice_roc(s,1),perm] = NaN;
                end
                
                if length(Zresp_left_sac(:,s)) >3 && length(Zresp_right_sac(:,s))>3
                    [sac_roc(s,1),~,perm] = rocN(Zresp_left_sac(:,s),Zresp_right_sac(:,s),100,1000);  %set the left as x axis
%                     [~,sac_roc(s,2)] = ttest2(Zresp_left_sac(:,s),Zresp_right_sac(:,s)); % Use ttest2 instead of permutation for speed (the outcomes are quite similar). HH20140601                    
                    sac_roc(s,2) = perm.pValue;
                    if ~isnan(perm.pValue)
                        shuffle_sac_roc(s,:) = perm.auROCPerm;
                    else
                        shuffle_sac_roc(s,:) = NaN(1000,1);
                    end
                else
                    [sac_roc(s,1),perm] = NaN;
                end     
                if s==7
                    shuffle_index = perm.randPermSortIndex;      % when trial end early, NaN value affect rocN use,  we use first time point returned value
                end
            end

            
            
            % select neuron according to 500~1500ms stimulus period
            if ialign==1
                if length(allZresp_left_choice) >3 && length(allZresp_right_choice)>3
                    [all_choice_roc(1,1),~,perm] = rocN(allZresp_left_choice,allZresp_right_choice,100,1000);  %set the left as x axis
                    all_choice_roc(1,2) = perm.pValue;
%                     [~,all_choice_roc(1,2)] = ttest2(allZresp_left_choice,allZresp_right_choice); % Use ttest2 instead of permutation for speed (the outcomes are quite similar). HH20140601
                    if ~isnan(all_choice_roc(1,2))
                       all_shuffle_choice_roc = perm.auROCPerm;
                    else
                        all_shuffle_choice_roc = NaN(1000,1);
                    end
                    permuteN = length(all_shuffle_choice_roc);
                    X = all_shuffle_choice_roc;
                    p = sum(X' > X, 2) / permuteN; % one tail
                    p(X < 0.5) = 1 - p(X < 0.5); % two tail
                    p = p * 2;
                    p(p > 1) = 1;
                    all_shuffle_choice_p = p';
                else
                    [all_choice_roc(1,1),perm] = NaN;
                end
                
                if length(allZresp_left_sac) >3 && length(allZresp_right_sac)>3
                    [all_sac_roc(1,1),~,perm] = rocN(allZresp_left_sac,allZresp_right_sac,100,1000);  %set the left as x axis
                    all_sac_roc(1,2) = perm.pValue;
%                     [~,all_sac_roc(1,2)] = ttest2(allZresp_left_sac,allZresp_right_sac); % Use ttest2 instead of permutation for speed (the outcomes are quite similar). HH20140601
                    if ~isnan(all_sac_roc(1,2))
                        all_shuffle_sac_roc = perm.auROCPerm;
                    else
                        all_shuffle_sac_roc = NaN(1000,1);
                    end
                    permuteN = length(all_shuffle_sac_roc);
                    X = all_shuffle_sac_roc;
                    p = sum(X' > X, 2) / permuteN; % one tail
                    p(X < 0.5) = 1 - p(X < 0.5); % two tail
                    p = p * 2;
                    p(p > 1) = 1;
                    all_shuffle_sac_p = p';
                else
                    [all_sac_roc(1,1),perm] = NaN;
                end
            end
            
            arranged_choice_mat = repmat(arranged_choice,1,1000);
            arranged_sac_mat = repmat(arranged_sac,1,1000);
            shuffle_choice = arranged_choice_mat(shuffle_index);
            shuffle_sac = arranged_sac_mat(shuffle_index);
%             spike_rates_mat = repmat(spike_rates{f,1}',1,1000);
%             z_fir_mat  = repmat(z_fir_per_trial,1,1000);

            
            [rr,pp] = corrcoef(relativeheading_per_trial{f,1}, spike_rates{f,1}');
            line(1,1) = rr(1,2);line(1,2) = pp(1,2);
            
            [rchoice,pchoice] = corrcoef(choice_per_trial{f,1}, z_fir_per_trial);
            choice_line(1,1) = rchoice(1,2);choice_line(1,2) = pchoice(1,2);
            
%             [shuffle_rchoice,shuffle_pchoice] = corrcoef(shuffle_seechoice, z_fir_per_trial);
%             shuffle_linechoice(1,1) = rchoice(1,2);shuffle_linechoice(1,2) = pchoice(1,2);
            arranged_choice_fir = [z_fir_per_trial(see_choice==-1,:);z_fir_per_trial(see_choice==1)];
            %check
            arranged_choice_fir_mat = repmat(arranged_choice_fir,1,1000);
            check = arranged_choice_fir_mat(shuffle_index);
            [shuffle_rchoice,shuffle_pchoice] = corr(check,arranged_choice);
%             [shuffle_rchoice,shuffle_pchoice] = corr(shuffle_choice,arranged_choice_fir); this is wrong
            shuffle_linechoice(:,1) = shuffle_rchoice; shuffle_linechoice(:,2) = shuffle_pchoice;
            

            [rsac,psac] = corrcoef(see_sac,spike_rates{f,1}');
            linesac(1,1) = rsac(1,2);linesac(1,2) = psac(1,2);
            
%             sac_fir = spike_rates{f,1}';
            sac_fir = fir_per_trial;
            arranged_sac_fir = [sac_fir(see_sac==-1,:);sac_fir(see_sac==1)];
            %check
            arranged_sac_fir_mat = repmat(arranged_sac_fir,1,1000);
            check = arranged_sac_fir_mat(shuffle_index);
            [shuffle_rsac,shuffle_psac] = corr(check,arranged_sac);
%             [shuffle_rsac,shuffle_psac] = corr(shuffle_sac,sac_fir); this is wrong
            shuffle_linesac(:,1) = shuffle_rsac; shuffle_linesac(:,2) = shuffle_psac;       
            
%             if  logical(line(1,1)>0)
%                 CP_T1L{f,ialign}(:,1) = 1-choice_roc_T1L(:,1);
%                 CP_T1R{f,ialign}(:,1) = 1-choice_roc_T1R(:,1);
%                 CP{f,ialign}(:,1) = 1-choice_roc(:,1);
%                 allCP{f,ialign}(:,1) = 1-all_choice_roc(:,1);
%             else
%                 CP_T1L{f,ialign}(:,1) = choice_roc_T1L(:,1);
%                 CP_T1R{f,ialign}(:,1) = choice_roc_T1R(:,1);
%                 CP{f,ialign}(:,1) = choice_roc(:,1);
%                 allCP{f,ialign}(:,1) = all_choice_roc(:,1);
%             end
%             CP_T1L{f,ialign}(:,2) = choice_roc_T1L(:,2);
%             CP_T1R{f,ialign}(:,2) = choice_roc_T1R(:,2);
%             CP{f,ialign}(:,2) = choice_roc(:,2);
%             allCP{f,ialign}(:,2) = all_choice_roc(:,2);
            if  logical(choice_line(1,1)>0)
%             if all_choice_roc(:,1)<0.5
                ROC_choice{f,ialign}(:,1) = 1-choice_roc(:,1);
                allROC_choice{f,ialign}(:,1) = 1-all_choice_roc(:,1);   %(all time period mean)
            else

                ROC_choice{f,ialign}(:,1) = choice_roc(:,1);
                allROC_choice{f,ialign}(:,1) = all_choice_roc(:,1);
            end
            ROC_choice{f,ialign}(:,2) = choice_roc(:,2);
            allROC_choice{f,ialign}(:,2) = all_choice_roc(:,2);
            

            shuffle_ROC_choice{f,ialign} = shuffle_choice_roc;
            shuffle_allROC_choice{f,ialign} = all_shuffle_choice_roc;
            choice_rightprefer = shuffle_linechoice(:,1)>0;
%             choice_rightprefer = shuffle_allROC_choice{f,ialign}<0.5;
            shuffle_ROC_choice{f,ialign}(:,choice_rightprefer) = 1-shuffle_choice_roc(:,choice_rightprefer);% time x shuffle_num
            shuffle_allROC_choice{f,ialign}(choice_rightprefer)  = 1-all_shuffle_choice_roc(choice_rightprefer);
            
            p_shuffle_allROC_choice{f,ialign} = all_shuffle_choice_p;%(1 x shufflenum) 
 
%             if all_sac_roc(:,1)<0.5
            if logical(linesac(1,1)>0) %prefer saccade rightward
                ROC_sac{f,ialign}(:,1) = 1-sac_roc(:,1);
                allROC_sac{f,ialign}(:,1) =1- all_sac_roc(:,1);
            else
                ROC_sac{f,ialign}(:,1) = sac_roc(:,1);
                allROC_sac{f,ialign}(:,1) = all_sac_roc(:,1);
            end
            ROC_sac{f,ialign}(:,2) = sac_roc(:,2);
            allROC_sac{f,ialign}(:,2) = all_sac_roc(:,2);
%             if logical(shuffle_linesac(1,1)>0) %prefer saccade rightward
%                 shuffle_ROC_sac{f,ialign}(:,1) = 1-shuffle_sac_roc(:,1);
%             else
%                 shuffle_ROC_sac{f,ialign}(:,1) = shuffle_sac_roc(:,1);
%             end
%             shuffle_ROC_sac{f,ialign}(:,2) = shuffle_sac_roc(:,2);
            shuffle_ROC_sac{f,ialign} = shuffle_sac_roc;
            shuffle_allROC_sac{f,ialign} = all_shuffle_sac_roc;
            sac_rightprefer = shuffle_linesac(:,1)>0;
%             sac_rightprefer = shuffle_allROC_sac{f,ialign}<0.5;
            shuffle_ROC_sac{f,ialign}(:,sac_rightprefer) = 1-shuffle_sac_roc(:,sac_rightprefer);% time x shuffle_num
            shuffle_allROC_sac{f,ialign}(sac_rightprefer)  = 1-all_shuffle_sac_roc(sac_rightprefer);
            
            p_shuffle_allROC_sac{f,ialign} = all_shuffle_sac_p;%(1 x shufflenum) 
            
            

            mat_ROC_choice{iarea,ialign}(f,:) = ROC_choice{f,ialign}(:,1);           
            mat_P_allROC_choice{iarea,ialign}(f,:) = allROC_choice{f,ialign}(:,2);
            for ishuffle = 1: 1000
                mat_shuffleROC_choice{iarea,ialign,ishuffle}(f,:) = shuffle_ROC_choice{f,ialign}(:,ishuffle);
            end

            mat_P_shuffleallROC_choice{iarea,ialign}(f,:) = p_shuffle_allROC_choice{f,ialign}; 

            mat_ROC_sac{iarea,ialign}(f,:) = ROC_sac{f,ialign}(:,1);
            mat_P_allROC_sac{iarea,ialign}(f,:) = allROC_sac{f,ialign}(:,2);
            for ishuffle = 1:1000
                mat_shuffleROC_sac{iarea,ialign,ishuffle}(f,:) = shuffle_ROC_sac{f,ialign}(:,ishuffle);
            end
            mat_P_shuffleallROC_sac{iarea,ialign}(f,:) = p_shuffle_allROC_sac{f,ialign}; 

        end
    select_allROC_choice{iarea,ialign} = mat_P_allROC_choice{iarea,1}<0.05;
    select_shuffleallROC_choice{iarea,ialign} =   mat_P_shuffleallROC_choice{iarea,1}<0.05;
    if iarea == 2
        number_choice_STS = sum(select_allROC_choice{1,1}(:));
        number_choice_LIP = sum(select_allROC_choice{2,1}(:));
    end
    select_allROC_sac{iarea,ialign} =   mat_P_allROC_sac{iarea,1}(:)<0.05;
    select_shuffleallROC_sac{iarea,ialign} =   mat_P_shuffleallROC_sac{iarea,1}<0.05;    
    if iarea==2
    number_saccade_STS = sum(select_allROC_sac{1,1});
    number_saccade_LIP = sum(select_allROC_sac{2,1});
    end
    mean_ROC_choice{iarea,ialign} = squeeze(nanmean(mat_ROC_choice{iarea,ialign}( select_allROC_choice{iarea,ialign},:)));
    fprintf('iarea=%d, ialign=%d, CP_num=%d\n', iarea, ialign, sum(select_allROC_choice{iarea,ialign}));

    sem_ROC_choice{iarea,ialign} = nanstd(mat_ROC_choice{iarea,ialign}(select_allROC_choice{iarea,ialign},:))/sqrt(size(mat_ROC_choice{iarea,ialign}( select_allROC_choice{iarea,ialign},:),1)); CI_ROC_choice{iarea,ialign} = sem_ROC_choice{iarea,ialign}*1.96;

%     tmp = cat(3, mat_shuffleROC_choice{iarea,ialign,:});
%     mean_shuffle_ROC_choice{iarea,ialign} = squeeze(nanmean(tmp, 1));
%     mmean_shuffle_ROC_choice{iarea,ialign} = nanmean(mean_shuffle_ROC_choice{iarea,ialign},2);
%     SD_shuffle_ROC_choice{iarea,ialign} = nanstd(mean_shuffle_ROC_choice{iarea,ialign},0,2); 
%     p_ROC_choice{iarea,ialign} = mean(mean_ROC_choice{iarea,ialign}' <= mean_shuffle_ROC_choice{iarea,ialign} ,2);
    
    for ishuffle = 1:1000
        mean_shuffle_ROC_choice{iarea,ialign}(ishuffle,:) = nanmean(mat_shuffleROC_choice{iarea,ialign,ishuffle}(select_shuffleallROC_choice{iarea,ialign}(:,ishuffle),:),1);
    end
    mmean_shuffle_ROC_choice{iarea,ialign} = nanmean(mean_shuffle_ROC_choice{iarea,ialign});
    SD_shuffle_ROC_choice{iarea,ialign} = nanstd(mean_shuffle_ROC_choice{iarea,ialign},0,1);
    CI_shuffle_ROC_choice{iarea,ialign} =  SD_shuffle_ROC_choice{iarea,ialign} ./1000*1.96;
    p_ROC_choice{iarea,ialign} = mean(mean_ROC_choice{iarea,ialign}<= mean_shuffle_ROC_choice{iarea,ialign} ,1);   
    
 
%     tmp = cat(3, mat_shuffleROC_sac{iarea,ialign,:});
%     mean_shuffle_ROC_sac{iarea,ialign} = squeeze(nanmean(tmp, 1));
%     mmean_shuffle_ROC_sac{iarea,ialign} = nanmean(mean_shuffle_ROC_sac{iarea,ialign},2);
%     SD_shuffle_ROC_sac{iarea,ialign} = nanstd(mean_shuffle_ROC_sac{iarea,ialign},0,2); 
%     p_ROC_sac{iarea,ialign} = mean(mean_ROC_sac{iarea,ialign}' <= mean_shuffle_ROC_sac{iarea,ialign} ,2);
    
 

    mean_ROC_sac{iarea,ialign} = squeeze(nanmean(mat_ROC_sac{iarea,ialign}( select_allROC_sac{iarea,ialign},:)));
    sem_ROC_sac{iarea,ialign} = nanstd(mat_ROC_sac{iarea,ialign}(select_allROC_sac{iarea,ialign},:))/sqrt(size(mat_ROC_sac{iarea,ialign}( select_allROC_sac{iarea,ialign},:),1));  CI_ROC_sac{iarea,ialign} = sem_ROC_sac{iarea,ialign}*1.96;
%     [~,P_ROCsac{iarea,ialign}] = ttest2(mat_ROCsac{iarea,ialign},mat_shuffleROCsac{iarea,ialign},'Tail','right');
    for ishuffle = 1:1000
        mean_shuffle_ROC_sac{iarea,ialign}(ishuffle,:) = nanmean(mat_shuffleROC_sac{iarea,ialign,ishuffle}(select_shuffleallROC_sac{iarea,ialign}(:,ishuffle),:),1);
    end
    mmean_shuffle_ROC_sac{iarea,ialign} = nanmean(mean_shuffle_ROC_sac{iarea,ialign});
    SD_shuffle_ROC_sac{iarea,ialign} = nanstd(mean_shuffle_ROC_sac{iarea,ialign},0,1);
    CI_shuffle_ROC_sac{iarea,ialign}  = SD_shuffle_ROC_sac{iarea,ialign}./1000*1.96;
    p_ROC_sac{iarea,ialign} = mean(mean_ROC_sac{iarea,ialign}<= mean_shuffle_ROC_sac{iarea,ialign} ,1);
    
    
    if iarea==2
        [~,P_area1vs2{1,ialign}] = ttest2(mat_ROC_choice{1,ialign},mat_ROC_choice{2,ialign});
        [~,Psac_area1vs2{1,ialign}] = ttest2(mat_ROC_sac{1,ialign},mat_ROC_sac{2,ialign});
    end
    end
    
%%
% lineColor = [58/255,98/255,205/255;238/255,48/255,167/255;];
if iarea==1
    saccolor=[30,114,255]/255; %saccade
else
    saccolor=[30,64,139]/255; %saccade
end
abschoicecolor=[139,10,80]/255; %choice
p_critical=0.05;
lineStyle = '-';
transparent = 1;

if iarea ==2
    figure(1)
    width = 100; height = 0.001;
    h_subplot1 = tight_subplot(1,2,[0.1 0.02],[0.1 0.1],[0.1 0.1],[],temp_duration_ratio);
for iiarea = 1:2
    lineColor = min([abschoicecolor(1,:)+[1,1,1]*1*(1-1*iiarea/2);1,1,1]);
    for ialign = 1:2
        set(gcf,'CurrentAxes',h_subplot1(ialign)); ax(ialign) = gca;
        h = shadedErrorBar(t_center_step{1,ialign},mean_ROC_choice{iiarea,ialign},CI_ROC_choice{iiarea,ialign},'lineprops',{'Color',lineColor,'LineStyle',lineStyle,'LineWidth',1.5},'transparent',transparent);
        h = shadedErrorBar(t_center_step{1,ialign},mmean_shuffle_ROC_choice{iiarea,ialign},SD_shuffle_ROC_choice{iiarea,ialign},'lineprops',{'Color',lineColor,'LineStyle','--','LineWidth',1.5},'transparent',transparent);
        set(h.mainLine,'LineWidth',2);  hold on;
        h_legend(iiarea) = h.mainLine;
        hold on;
        plotsignificant(t_center_step{1,ialign}(p_ROC_choice{iiarea,ialign}<p_critical),0.46*ones(size(t_center_step{1,ialign}(p_ROC_choice{iiarea,ialign}<p_critical)))+iiarea*0.004,width, height,lineColor)
%         plot(t_center_step{1,ialign}(P_CP{iiarea,ialign}<p_critical),(0.45*ones(size(t_center_step{1,ialign}(P_CP{iiarea,ialign}<p_critical))))+iiarea*0.01,'s','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',5)
        if iiarea==2
            
            plotsignificant(t_center_step{1,ialign}(P_area1vs2{1,ialign}<p_critical),0.46*ones(size(t_center_step{1,ialign}(P_area1vs2{1,ialign}<p_critical))),width, height,[0.5 0.5 0.5])            
%             plot(t_center_step{1,ialign}(P_area1vs2{1,ialign}<p_critical),0.45*ones(size(t_center_step{1,ialign}(P_area1vs2{1,ialign}<p_critical))),'s','Color','k','MarkerFaceColor','k','MarkerSize',5)

        %     linkaxes(ax,'y')
        
            set(gcf,'CurrentAxes',h_subplot1(ialign))
            if ialign==1
                xlim([-200 1600])
                ylim([0.45,0.58])
                stim_on = 0; stim_off = align_offsets_others{1,1}{1,1}(1,1);
                plot([stim_on,stim_on],ylim,'k','linewidth',2)
                plot([stim_off,stim_off],ylim,'k','linewidth',2)
                xticks([0,750,1500])
                yticks([0.45,0.5,0.55])
                xticklabels({'stim on', '750','stim off'})
                title('CP(abstract choice)')
                if iiarea==2
                    txt_legend = {'IPS','STS'};
    %                 legend(ax(1),h_legend,txt_legend,'Location','best','color','none','AutoUpdate','off');
                end
            else
                xlim([-300,200])
                ylim([0.45,0.58])
                sac_on = 0;
                stim_off = [];
                for ifile = 1:length(align_offsets_others)
                    stim_off =[ stim_off;align_offsets_others{ifile,1}{1,2}(:,2)];
                end
                meanplusstd_stim_off = mean(stim_off)+std(stim_off);
                meanminusstd_stim_off = mean(stim_off)-std(stim_off);
                plot([sac_on,sac_on],ylim,'k','linewidth',2)
                plot([meanplusstd_stim_off,meanplusstd_stim_off],ylim,'k--','linewidth',1.5)
                plot([meanminusstd_stim_off,meanminusstd_stim_off],ylim,'k--','linewidth',1.5)
                xticks([0])
                xticklabels({'saccade on'})
                set(gca,'YColor','none')

            end
        end
    end
    SetFigure
end

figure(2)
width = 100; height = 0.001;
h_subplot2 = tight_subplot(1,2,[0.1 0.02],[0.1 0.1],[0.1 0.1],[],temp_duration_ratio);
for iiarea = 1:2
    if iiarea==1
        saccolor=[30,114,255]/255; %saccade
    else
        saccolor=[30,64,139]/255; %saccade
    end
%     lineColor = min([saccolor(1,:)+[1,1,1]*(1-1*iiarea/2);1,1,1]);
    lineColor = saccolor;
    for ialign = 1:2
        set(gcf,'CurrentAxes',h_subplot2(ialign)); ax(ialign) = gca;
        h = shadedErrorBar(t_center_step{1,ialign},mean_ROC_sac{iiarea,ialign},CI_ROC_sac{iiarea,ialign},'lineprops',{'Color',lineColor,'LineStyle',lineStyle','LineWidth',1.5},'transparent',transparent);
        h = shadedErrorBar(t_center_step{1,ialign},mmean_shuffle_ROC_sac{iiarea,ialign},SD_shuffle_ROC_sac{iiarea,ialign},'lineprops',{'Color',lineColor,'LineStyle','--','LineWidth',1.5},'transparent',transparent);
        set(h.mainLine,'LineWidth',2);  hold on;
        h_legend(iiarea) = h.mainLine;
        hold on;
        plotsignificant(t_center_step{1,ialign}(p_ROC_sac{iiarea,ialign}<p_critical),0.48*ones(size(t_center_step{1,ialign}(p_ROC_sac{iiarea,ialign}<p_critical)))+iiarea*0.005,width, height,lineColor)
%         plot(t_center_step{1,ialign}(P_ROCsac{iiarea,ialign}<p_critical),(0.45*ones(size(t_center_step{1,ialign}(P_ROCsac{iiarea,ialign}<p_critical))))+iiarea*0.01,'s','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',5)
        if iiarea==2      
            hold on;
            plotsignificant(t_center_step{1,ialign}(Psac_area1vs2{1,ialign}<p_critical),0.48*ones(size(t_center_step{1,ialign}(Psac_area1vs2{1,ialign}<p_critical))),width, height, [0.5,0.5,0.5])
%             plot(t_center_step{1,ialign}(Psac_area1vs2{1,ialign}<p_critical),0.45*ones(size(t_center_step{1,ialign}(Psac_area1vs2{1,ialign}<p_critical))),'s','Color','k','MarkerFaceColor','k','MarkerSize',5)      

        %     linkaxes(ax,'y')
        
            set(gcf,'CurrentAxes',h_subplot2(ialign))
            if ialign==1
                xlim([-200 1600])
                ylim([0.45,0.68])
                stim_on = 0; stim_off = align_offsets_others{1,1}{1,1}(1,1);
                plot([stim_on,stim_on],ylim,'k','linewidth',2)
                plot([stim_off,stim_off],ylim,'k','linewidth',2)
                xticks([0,750,1500])
                yticks([0.45,0.55,0.65])
    %             xticklabels({'stim on', '500','1000','stim off'})
                title('Saccade ROC')
                if iiarea==2
                    txt_legend = {'IPS','STS'};
    %                 legend(ax(1),h_legend,txt_legend,'Location','best','color','none','AutoUpdate','off');
                end
            else
                xlim([-300,200])
                ylim([0.45,0.68])
                sac_on = 0;
                stim_off = [];
                for ifile = 1:length(align_offsets_others)
                    stim_off =[ stim_off;align_offsets_others{ifile,1}{1,2}(:,2)];
                end
                meanplusstd_stim_off = mean(stim_off)+std(stim_off);
                meanminusstd_stim_off = mean(stim_off)-std(stim_off);
                plot([sac_on,sac_on],ylim,'k','linewidth',2)
                plot([meanplusstd_stim_off,meanplusstd_stim_off],ylim,'k--','linewidth',1.5)
                plot([meanminusstd_stim_off,meanminusstd_stim_off],ylim,'k--','linewidth',1.5)
                xticks([0])
                xticklabels({'saccade on'})
                set(gca,'YColor','none')
            end
       end
    end
    SetFigure
end
end
end