% rerange neural data for mTDR analysis 2024/10/27
% generate task conditions X, size(X) = (Nmax,var_unique)
% generate neural response Y, size(Y) = (Nneuron,time_len,Nmax)
function data_rerange(FolderPath_colorHD, time_bin, savepath, savename, shuffle_trial,randseed)
if nargin==0
    %% load all cell color task data
    % FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Fridaycolor113up_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\Fridaycolor113up_NP\',...
    %     'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Fridaycolor113down_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\Fridaycolor113down_NP\'};
    % FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Fridaycolor113down_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\Fridaycolor113down_NP\'};
    % FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorNTSTS_LA\'};
    % FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorSTS_NP\'};
    % FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorSTS_LA\'};
    % FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseSTS_LA\';};
    % FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpost\DrogoncolorcoarseLIP\';'Z:\Data\Tempo\Batch\fingerpost\DrogoncolorcoarseLIPsimultaneous\';'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseLIP_LA\'};
    % FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorSTS_LA\';'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorSTS_NP\'};
    FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseIPS_LA\';};
    time_bin = 100;%ms
    savename = ['F_IPScolor_shuffle_all_',num2str(time_bin),'ms'];
    shuffle_trial = 0;

end


nheading = 5; 
maxtrial_num = 60;

if length(FolderPath_colorHD)==1   
    Files_colorHD = dir(fullfile(FolderPath_colorHD{1},'*PSTH.mat'));
    pathfilenum_color = length(Files_colorHD);
else
    Files_colorHD = [];
    for ipath = 1:length(FolderPath_colorHD)
        ipathfiles = dir(fullfile(FolderPath_colorHD{ipath},'*PSTH.mat'));
        Files_colorHD = [Files_colorHD;ipathfiles];
        pathfilenum_color(ipath) = length(ipathfiles);
    end
    pathfilenum_color = cumsum(pathfilenum_color);
end
File_work_count = 0;
for ct = 1:length(Files_colorHD)
    File_work_count = File_work_count+1;
    current_path_pod = find(pathfilenum_color<File_work_count);
    if isempty(current_path_pod)
        current_path = 1;
    else
        current_path = max(current_path_pod) + 1;
    end
    list_color{File_work_count} = load([FolderPath_colorHD{current_path} Files_colorHD(ct).name]);
    cell_channel = list_color{1,ct}.result.SpikeChan;
    TF = isstrprop(list_color{1,ct}.result.FILE,'digit');
    charloc = find(TF==0);
    monkeyID = str2num(list_color{1,ct}.result.FILE(charloc(1)+1:charloc(2)-1));
    cellnum_color(ct) = str2double(list_color{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;
    blocknum_color(ct) = str2double(list_color{1,ct}.result.FILE(charloc(3)+1:end));
    cellID_color(ct) = str2double(list_color{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;   %+blocknum(ct)/10;
end
%% load all cell non color HD task data(only for get neural recorded in both task)
% FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridayMST_NP\'};
% FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogonIPS_LA\'};
FolderPath_HD = [];
if ~isempty(FolderPath_HD )
    if length(FolderPath_HD)==1
        Files_HD = dir(fullfile(FolderPath_HD{1},'*PSTH.mat'));
        pathfilenum_HD = length(Files_HD);
    else
        Files_HD = [];
        for ipath = 1:length(FolderPath_HD)
            ipathfiles = dir(fullfile(FolderPath_HD{ipath},'*PSTH.mat'));
            Files_HD = [Files_HD;ipathfiles];
            pathfilenum_HD(ipath) = length(ipathfiles);
        end
        pathfilenum_HD = cumsum(pathfilenum_HD);
    end
    File_work_count = 0;
    for ct = 1:length(Files_HD)
        File_work_count = File_work_count+1;
        current_path_pod = find(pathfilenum_HD<File_work_count);
        if isempty(current_path_pod)
            current_path = 1;
        else
            current_path = max(current_path_pod) + 1;
        end
        list_HD{File_work_count} = load([FolderPath_HD{current_path} Files_HD(ct).name]);
        cell_channel = list_HD{1,ct}.result.SpikeChan;
        TF = isstrprop(list_HD{1,ct}.result.FILE,'digit');
        charloc = find(TF==0);
        monkeyID = str2num(list_HD{1,ct}.result.FILE(charloc(1)+1:charloc(2)-1));
        cellnum_HD(ct) = str2double(list_HD{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;
        blocknum_HD(ct) = str2double(list_HD{1,ct}.result.FILE(charloc(3)+1:end));
        cellID_HD(ct) = str2double(list_HD{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;   %+blocknum(ct)/10;
    end
end
%%  select specific cells(have significant 1D tuning or mem sac tuning)
%cellID_1D
%     FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dMST_LA\','Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dMST_LA\'};
% FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dVIP_LA\','Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dLIP_LA\';};
% FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dSTS_LA\','Z:\Data\Tempo\Batch\fingerpostFriday_NP\Friday1dNTSTS_NP\'};
%     FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dIPS_LA\'};
FolderPath_1D = [];
if ~isempty(FolderPath_1D )
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
end
%%
%     %cellID_mem
%    FolderPath_mem = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridayMemIPS_LA\'};
% FolderPath_mem = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaymemSTS_LA\'};
FolderPath_mem = [];
if ~isempty(FolderPath_mem)
    if length(FolderPath_mem)==1
        Files_mem = dir(fullfile(FolderPath_mem{1},'*MemSac.mat'));
        pathfilenum = length(Files_mem);
    else
        Files_mem = [];
        pathfilenum = [];
        for ipath = 1:length(FolderPath_mem)
            ipathfiles = dir(fullfile(FolderPath_mem{ipath},'*MemSac.mat'));
            Files_mem = [Files_mem;ipathfiles];
            pathfilenum(ipath) = length(ipathfiles);
        end
        pathfilenum = cumsum(pathfilenum);
    end

    File_work_count_mem = 0;
    for ct = 1:length(Files_mem)
        File_work_count_mem = File_work_count_mem+1;
        current_path_pod = find(pathfilenum<File_work_count_mem);
        if isempty(current_path_pod)
            current_path = 1;
        else
            current_path = max(current_path_pod) + 1;
        end
        list_mem{File_work_count_mem} = load([FolderPath_mem{current_path} Files_mem(ct).name]);
        p_mem(ct) = min(list_mem{1,ct}.result.p);
        cell_channel = list_mem{1,ct}.result.SpikeChan;
        TF = isstrprop(list_mem{1,ct}.result.FILE,'digit');
        charloc = find(TF==0);
        monkeyID = str2num(list_mem{1,ct}.result.FILE(charloc(1)+1:charloc(2)-1));
        cellnum_mem(ct) = str2double(list_mem{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;
        blocknum_mem(ct) = str2double(list_mem{1,ct}.result.FILE(charloc(3)+1:end));
        cellID_mem(ct) = str2double(list_mem{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;   %+blocknum(ct)/10;
    end
    select_cellID_mem = cellID_mem(p_mem<0.05);
end
%%
% [selective_one,~,~] = union(select_cellID_1D,select_cellID_mem);
% [cellID_both,cellindex_color,cellindex_HD] = intersect(cellID_color,cellID_HD); % select cell record  for both task
% [cellID_pick,cellindex,cellindex_1D] = intersect(cellID_both,select_cellID_1D);

%     [cellID_pick,cellindex,cellindex_1D] = intersect(cellID_color,select_cellID_1D);
[cellID_pick,cellindex,cellindex] = intersect(cellID_color,cellID_color);
%%
icount = 0; % count the concatenated neuron num
Ymat = [];
% list = list_color(cellindex_color);
list = list_color;
% cellID_pick = cellID_both;
% cellindex = cellindex_color;
for ct = 1:length(cellID_pick)
    cell_channel = list{1,cellindex(ct)}.result.SpikeChan; 
    PSTH_cell{ct} = list{1,cellindex(ct)}.result.spike_hist(1,1);
    heading_per_cell{ct} = list{1,cellindex(ct)}.result.heading_per_trial;
    targetloc_per_cell{ct} = list{1,cellindex(ct)}.result.targetloc_per_trial;
    choice_per_cell{ct} = list{1,cellindex(ct)}.result.choice_per_trial;
    coherence_per_cell{ct} = list{1,cellindex(ct)}.result.coherence_per_trial;
    spike_aligned = list{1,cellindex(ct)}.result.spike_aligned;
%     t_centers_per_cell{ct} = list(1,cellindex(ct)).result.t_centers{1,1};
    timemarker_per_cell{ct} = list{1,cellindex(ct)}.result.align_offsets_others{1,1};
    correct_per_cell{ct} = sign(heading_per_cell{ct}).*sign(choice_per_cell{ct});
    sacdir_per_cell{ct} = targetloc_per_cell{ct}.*choice_per_cell{ct};
    align_markers = list{1,cellindex(ct)}.result.align_markers;
%     PSTH_per_trial = cell2mat(PSTH_cell{1,ct});
    coherence_per_trial = coherence_per_cell{ct};
    if length(unique(coherence_per_trial)) == 1
        heading_per_trial = heading_per_cell{ct};
    else
        heading_per_trial = sign(heading_per_cell{ct}).*coherence_per_trial;
    end
    targetloc_per_trial = sign(targetloc_per_cell{ct});
    choice_per_trial = choice_per_cell{ct};   
    sacdir_per_trial = targetloc_per_trial.*choice_per_trial; %saccde left mark -1, saccade right mark 1
    lastsacdir_per_trial = [0,sacdir_per_trial(1:end-1)];
    stimend_per_trial = timemarker_per_cell{ct}(:,1);
    correct_per_trial = correct_per_cell{ct};
    
    % calculate each neural prefer
    stim_on = find(list{1,cellindex(ct)}.result.spike_aligned{2,1}>0,1,'first');
    stim_off = find(list{1,cellindex(ct)}.result.spike_aligned{2,1}<1500,1,'last');
    firing_rate = sum(list{1,cellindex(ct)}.result.spike_aligned{1,1}(:,stim_on:stim_off),2)/1500;
    %heading prefer(need Zscore)
    z_choice_firing = [];
    unique_choice = unique(choice_per_trial);
    for ichoice = 1:length(unique_choice)
        select = choice_per_trial==unique_choice(ichoice);
        z_fir = firing_rate(select);
        z_fir = (z_fir-mean(z_fir))/std(z_fir);
        z_choice_firing(select) = z_fir;         
    end
    [r_heading,p_heading] = corrcoef(heading_per_trial,z_choice_firing);   
    if r_heading>0
        prefer_heading_sign = 1;
    else
        prefer_heading_sign = -1;
    end        
    
    %choice prefer(need Zscore)
    z_head_firing = [];
    unique_heading = unique(heading_per_trial);
    for iheading = 1:length(unique_heading)
        select = heading_per_trial==unique_heading(iheading);
        z_fir = firing_rate(select);
        z_fir = (z_fir-mean(z_fir))/std(z_fir);
        z_head_firing(select) = z_fir;    
    end
    [r_choice,p_choice] = corrcoef(choice_per_trial,z_head_firing);   
    if r_choice>0
        prefer_choice_sign = 1;
    else
        prefer_choice_sign = -1;
    end
        
    %target location prefer
    [r_targloc,p_targloc] = corrcoef(targetloc_per_trial,firing_rate);
    if r_targloc>0
        prefer_targloc_sign = 1;
    else
        prefer_targloc_sign = -1;
    end
    %saccade direction prefer
    [r_sacdir,p_sacdir] = corrcoef(sacdir_per_trial,firing_rate);
    if r_sacdir>0
        prefer_sacdir_sign = 1;
    else
        prefer_sacdir_sign = -1;
    end
    %last saccade direction prefer
    [r_lastsacdir,p_sacdir] = corrcoef(lastsacdir_per_trial,firing_rate);
    if r_lastsacdir>0
        prefer_lastsacdir_sign = 1;
    else
        prefer_lastsacdir_sign = -1;
    end
    %%Calculate PSTH using sliding windows
    % Anne Churchland NN, bin = 10 ms with Gaussian filter (sigma = 50 ms)
    if time_bin == 100
        binSize_rate = 100;  % in ms
        stepSize_rate = 100; % in ms
        smoothFactor = 125; % in ms !!
    elseif time_bin == 20
        binSize_rate = 20;  % in ms
        stepSize_rate = 20; % in ms
        smoothFactor = 100; % in ms !!
    end

    spike_timeWin = 1;

    for j = 1   % For each desired marker        
        if shuffle_trial == 1
            random_index_time = randperm(length(spike_aligned{1,j}));
            spike_aligned{1,j} = spike_aligned{1,j}(:,random_index_time);
        end
        t_centers{j} = align_markers{j,2} + binSize_rate/2 : stepSize_rate : align_markers{j,3} - binSize_rate/2; % Centers of PSTH time windows
        
        spike_hist{j} = zeros(size(spike_aligned{1,j},1),length(t_centers{j})); % Preallocation
        for k = 1:length(t_centers{j})
            winBeg = ceil(((k-1) * stepSize_rate) / spike_timeWin) + 1;
            winEnd = ceil(((k-1) * stepSize_rate + binSize_rate) / spike_timeWin) + 1;
            spike_hist{j}(:,k) = sum(spike_aligned{1,j}(: , winBeg:winEnd),2) / binSize_rate*1000 ;  % in Hz
        end
        
        % --- Smoothing ---
        % Anne Churchland NN, bin = 10 ms with Gaussian filter (sigma = 50 ms)
        % Note that this only influence PSTH calculation (spike_hist), not CP
        if smoothFactor > 0
            for i = 1:size(spike_hist{j},1) % Each trial
                spike_hist{j}(i,:) =  GaussSmooth(t_centers{j},spike_hist{j}(i,:),smoothFactor);
            end
        end
    end
    PSTH_per_trial = spike_hist{1};
   
    %% sort trial according to task variable 
    select_unique_targetloc = [-1,1];

    select_unique_choice = unique(choice_per_trial);
    
    unique_heading = sort(unique(heading_per_trial));
    if nheading==7
        relative_unique_heading = [-3,-2,-1,0,1,2,3]; % only for drogon some recording has 7 heading angle
    elseif nheading==5
        relative_unique_heading = [-2,-1,0,1,2]; % only for drogon some recording has 7 heading angle
    end
    
%     relative_unique_heading = [-4,-3,-2,-1,0,1,2,3,4];

    heading0_index = find(unique_heading == 0);
%     if isempty(heading0_index)
%        heading0_index = length(unique_heading)/2+1;
%        unique_heading = [unique_heading(1:length(unique_heading)/2),0,unique_heading(length(unique_heading)/2+1:end)];
%     end
    
    select_unique_heading = unique_heading;
    select_relative_heading = relative_unique_heading;
    select_unique_sac = [-1,1];  
    select_unique_lastsac = [-1,1];
    Xmat = [];
    each_Ymat = [];
    for iheading = 1:length(select_relative_heading)
        for itarget_loc = 1:length(select_unique_targetloc)
            for ichoice = 1:length(select_unique_choice)
                for ilastsac = 1:length(select_unique_choice)
                    curr_relative_heading = select_relative_heading(iheading);

                    if length(select_unique_heading)==7
                        if nheading == 7
                            curr_heading = select_unique_heading(heading0_index+curr_relative_heading);
                        elseif nheading == 5
                            if iheading == 1
                                curr_heading = [select_unique_heading(heading0_index+curr_relative_heading-1),select_unique_heading(heading0_index+curr_relative_heading)];
                            elseif ismember(iheading,[2,3,4])
                                curr_heading = select_unique_heading(heading0_index+curr_relative_heading);
                            elseif iheading == 5
                                curr_heading = [select_unique_heading(heading0_index+curr_relative_heading+1),select_unique_heading(heading0_index+curr_relative_heading)];
                            end
                        end
                    elseif length(select_unique_heading)==9
                        if nheading == 7
                            if iheading==1
                                curr_heading = [select_unique_heading(heading0_index+curr_relative_heading-1),select_unique_heading(heading0_index+curr_relative_heading)];
                            elseif ismember(iheading,[2,3,4,5,6])
                                curr_heading = select_unique_heading(heading0_index+curr_relative_heading);
                             elseif iheading ==7
                                curr_heading = [select_unique_heading(heading0_index+curr_relative_heading+1),select_unique_heading(heading0_index+curr_relative_heading)];
                            end
                        elseif nheading == 5
                            if iheading == 1
                                curr_heading = [select_unique_heading(heading0_index+curr_relative_heading-2),select_unique_heading(heading0_index+curr_relative_heading-1),select_unique_heading(heading0_index+curr_relative_heading)];
                            elseif ismember(iheading,[2,3,4])
                                curr_heading = select_unique_heading(heading0_index+curr_relative_heading);
                            elseif iheading == 5
                                curr_heading = [select_unique_heading(heading0_index+curr_relative_heading+2),select_unique_heading(heading0_index+curr_relative_heading+1),select_unique_heading(heading0_index+curr_relative_heading)];

                            end
                        end
                    end
                    curr_target_loc = select_unique_targetloc(itarget_loc);%*prefer_targloc_sign;
%                     curr_choice = select_unique_choice(ichoice)*prefer_choice_sign;
                    curr_choice = select_unique_choice(ichoice);%*prefer_heading_sign;                  
                    curr_sacdir = curr_target_loc*curr_choice;%*prefer_sacdir_sign;
                    last_sacdir = select_unique_lastsac(ilastsac);%*prefer_sacdir_sign;% use current saccade prefer direction as standard
                    curr_absheading = abs(curr_relative_heading);
                    rng(randseed)
                    random_index = randperm(length(choice_per_trial));
                    rand_choice_per_trial = choice_per_trial(random_index);
                    rand_heading_per_trial = heading_per_trial(random_index);
                    rand_targetloc_per_trial = targetloc_per_trial(random_index);
                    rand_lastsacdir_per_trial = lastsacdir_per_trial(random_index);
                    if shuffle_trial == 0
                        select = ismember(heading_per_trial,curr_heading)&targetloc_per_trial==curr_target_loc&...
                            choice_per_trial==curr_choice&lastsacdir_per_trial==last_sacdir;
                    elseif shuffle_trial == 1
                        select = ismember(rand_heading_per_trial,curr_heading)&rand_targetloc_per_trial==curr_target_loc&...
                            rand_choice_per_trial==curr_choice&rand_lastsacdir_per_trial==last_sacdir;                 
                    end
                    curr_relative_heading = sign(curr_relative_heading);
                    current_X = [curr_relative_heading,curr_target_loc,curr_choice,curr_sacdir,last_sacdir,curr_absheading,1];% the last P is time component
                    %                 current_X = [curr_target_loc,curr_choice,curr_sacdir,curr_absheading,1];% the last P is time component
                    current_Y = PSTH_per_trial(select,:);
                    current_Xmat = repmat(current_X,maxtrial_num,1);
                    [~,time_len] = size(PSTH_per_trial);
                    current_Ymat = zeros(time_len,maxtrial_num);
                    current_Ymat(:,1:sum(select)) =current_Y';
                    column_rank = randperm(size(current_Ymat,2));
                    current_Ymat  = current_Ymat(:,column_rank);
                    Xmat = [Xmat;current_Xmat];
                    each_Ymat = [each_Ymat,current_Ymat];
                    each_hk_sum = sum(each_Ymat);
                    each_hk = each_hk_sum~=0;
                end
            end
        end
    end
    Ymat(ct,:,:) = each_Ymat;
    hk(:,ct) = each_hk;
end

save([savepath savename],'Xmat','Ymat','t_centers','hk','cellID_pick','-v7.3')    
end
    