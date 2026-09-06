%% plot noise correlation for different color target location
2026/1/12
clear; close all;

% FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseSTS_LA\'};
% FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dSTS_LA\'};
% FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorSTS_LA\'};
% FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dSTS_LA\'};
% FolderPath = 'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTIPS_NP\';
% FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTMST_NP\'};
% FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostFriday_NP\Friday1dMST_NP\'};
FolderPath_HD = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseSTS_LA\';'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorSTS_LA\'};
FolderPath_1D = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dSTS_LA\';'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dSTS_LA\'};
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
    blocknum{ct} = list{1,ct}.result.FILE;
    cellID(ct) = str2double(list{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;   %+blocknum(ct)/10;
end
unique_session_name_color = unique(blocknum);
Files_session_color = cell(length(unique_session_name_color),1);
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
% [cellID_pick,cellindex,cellindex_1D] = intersect(cellID,select_cellID_1D);
[cellID_pick,cellindex,~] = intersect(cellID,cellID);
for f = 1:length(cellID_pick)
    neuronID_session = list{1,cellindex(f)}.result.FILE;
    neuronID_channel = list{1,cellindex(f)}.result.SpikeChan;
    neuronID =  [neuronID_session,'_',num2str(neuronID_channel),'_','*PSTH.mat'];
%     Files_color = dir(fullfile([FolderPath,neuronID_session,'r','*','_',neuronID_channel,'_','*PSTH.mat']));
   
    for isession = 1:length(unique_session_name_color)
        if strcmp(neuronID_session ,unique_session_name_color{isession})
            Files_session_color{isession,1}{end+1} = neuronID;
        end
    end
end

for isession = 1:length(unique_session_name_color)
if length(Files_session_color{isession,1}) < 2
    continue
end
list = [];tunning_per_cell_1 = [];tunning_per_cell_2 = []; fir_1 = [];fir_2 = [];tunning_per_cell=[]; fir_1_sac = [];
fir_1_nosac= []; fir_2_sac = []; fir_2_nosac= [];
    for f = 1:length(Files_session_color{isession,1})
        session_infor = Files_session_color{isession,1}{1,1};
        TF = isstrprop(session_infor,'digit');
        charloc = find(TF==0);
        monkeyID = str2double(session_infor(charloc(1)+1:charloc(2)-1));
        if monkeyID == 17
            pattern = fullfile(FolderPath_HD{1}, Files_session_color{isession}{f});
        elseif monkeyID == 24
            pattern = fullfile(FolderPath_HD{2}, Files_session_color{isession}{f});
        end
        fileinfo = dir(pattern);
        fullpath = fullfile(fileinfo(1).folder, fileinfo(1).name);
        list{f,1} = load(fullpath); % the list contain files recorded from the same session
    
        spike_hist = [{list{f,1}.result.spike_hist{1,1}},list{f,1}.result.t_centers{1,1}'];
        spike_aligned = [{list{f,1}.result.spike_aligned{1,1}},{list{f,1}.result.spike_aligned{2,1}}];

        heading_per_trial = (list{f,1}.result.heading_per_trial)';
        coherence_per_trial = (list{f,1}.result.coherence_per_trial)';
        if length(unique(coherence_per_trial)) >1
            relativeheading_per_trial = sign(heading_per_trial).*coherence_per_trial;
        else
            relativeheading_per_trial = heading_per_trial;
        end
        unique_relativeheading = unique(relativeheading_per_trial);
        choice_per_trial = (list{f,1}.result.choice_per_trial)';
        T1loc_per_trial = (list{f,1}.result.targetloc_per_trial)';
        saccade_per_trial = ((T1loc_per_trial)'.*list{f,1}.result.choice_per_trial)';
        saccade_per_trial = sign(saccade_per_trial);

            % check whether the T1 loc change will affect signal noise corr
        select_condition1 = sign(T1loc_per_trial)==-1; % select T1 on left
        select_condition2 = sign(T1loc_per_trial)== 1; % select T1 on right

%         select_condition1 = sign(saccade_per_trial)==-1; % select T1 on left
%         select_condition2 = sign(saccade_per_trial)== 1; % select T1 on right

        %for time window selection
        stim_dur = list{f,1}.result.align_offsets_others{1,2}(1,2) - list{1,1}.result.align_offsets_others{1,2}(1,1);
        stimulus_on= 0;
        stimulus_end = 0+stim_dur;


        %method1 calcalate with psth data        
        hist_time_window1(f,1) = find(spike_hist{1,2}(:) > stimulus_on+500,1);
        hist_time_window1(f,2) = find(spike_hist{1,2}(:) < stimulus_end+1500,1,'last');
        %method2 calculate with 0/1 spike data
        bin_time_window1(f,1) = find(spike_aligned{1,2}(:) > stimulus_on+0,1);
        bin_time_window1(f,2) = find(spike_aligned{1,2}(:) < stimulus_end+1500,1,'last');

        spike_per_cell = mean(spike_hist{1,1}(:,hist_time_window1(f,1):hist_time_window1(f,2)),2);
        
        % roc calculate saccade choice selectivity
        resp_left_sac =  spike_per_cell(saccade_per_trial ==-1,:);
        resp_right_sac = spike_per_cell(saccade_per_trial ==1,:);
        if length(resp_left_sac) >3 && length(resp_right_sac)>3
            [sac_roc,~,perm] = rocN(resp_left_sac,resp_right_sac,100,1000);  %set the left as x axis
            pvalue_sac_roc = perm.pValue; 
    %        [~,pvalue_sac_roc(f,1)] = ttest2(Zresp_left_sac,Zresp_right_sac);
        else
            [sac_roc,~,perm] = NaN;
        end

        % ------------- zscore for per cell----------------
        % zscore for heading
        fir_1{1,f} = [];
        fir_2{1,f} = [];

        unique_heading = unique(relativeheading_per_trial);
        for iheading = 1:length(unique_heading)
            select = (relativeheading_per_trial==unique_heading(iheading));
            select1 = (relativeheading_per_trial==unique_heading(iheading))&select_condition1;
            select2 = (relativeheading_per_trial==unique_heading(iheading))&select_condition2;
            mean_fir = mean(spike_per_cell(select));
%             mean_fir_1 = mean(spike_per_cell(select1));
%             mean_fir_2 = mean(spike_per_cell(select2));
%             tunning_per_cell_1{1,f}(iheading,1) = mean_fir_1';
%             tunning_per_cell_2{1,f}(iheading,1) = mean_fir_2';
            tunning_per_cell{1,f}(iheading,1) = mean_fir';

            fir_1{1,f} = [fir_1{1,f};spike_per_cell(select1)-mean_fir];
            fir_2{1,f} = [fir_2{1,f};spike_per_cell(select2)-mean_fir];
        end
            if pvalue_sac_roc<0.05
                fir_1_sac = [fir_1_sac,fir_1(1,f)];
                fir_2_sac = [fir_2_sac,fir_2(1,f)];
            else
                fir_1_nosac = [fir_1_nosac,fir_1(1,f)];
                fir_2_nosac = [fir_2_nosac,fir_2(1,f)];
            end
    end
%     units_fir_mat_1 = (cell2mat(tunning_per_cell_1));
%     sig_corr_matrix_1{isession} = corr(units_fir_mat_1,units_fir_mat_1);
%     units_fir_mat_2 = (cell2mat(tunning_per_cell_2));
%     sig_corr_matrix_2{isession} = corr(units_fir_mat_2,units_fir_mat_2);
    units_fir_mat = (cell2mat(tunning_per_cell));
    sig_corr_matrix{isession} = corr(units_fir_mat,units_fir_mat);
    
    fir_1_mat = cell2mat(fir_1);
    fir_2_mat = cell2mat(fir_2);
    noise_corr_matrix1{isession} = corr(fir_1_mat,fir_1_mat);
    noise_corr_matrix2{isession} = corr(fir_2_mat,fir_2_mat);
    if ~isempty(fir_1_sac)&&~isempty(fir_1_nosac)
        fir_1_sacmat = cell2mat(fir_1_sac);
        fir_2_sacmat = cell2mat(fir_2_sac);
        fir_1_nosacmat = cell2mat(fir_1_nosac);
        fir_2_nosacmat = cell2mat(fir_2_nosac);
        noise_corr_matrix1_sac{isession} = corr(fir_1_sacmat,fir_1_nosacmat);
        noise_corr_matrix2_sac{isession} = corr(fir_2_sacmat,fir_2_nosacmat);
        
        noise_corr_all_1_sac{isession,1} = tril(noise_corr_matrix1_sac{isession},-1);
        noise_corr_all_1_sac{isession,1} = noise_corr_all_1_sac{isession,1}(find(noise_corr_all_1_sac{isession,1}));
        noise_corr_all_2_sac{isession,1} = tril(noise_corr_matrix2_sac{isession},-1);
        noise_corr_all_2_sac{isession,1} = noise_corr_all_2_sac{isession,1}(find(noise_corr_all_2_sac{isession,1}));
    end
    
%     sig_corr_all_1{isession,1} = tril(sig_corr_matrix_1{isession},-1);
%     sig_corr_all_1{isession,1} = sig_corr_all_1{isession,1}(find(sig_corr_all_1{isession,1}));
%     sig_corr_all_2{isession,1} = tril(sig_corr_matrix_2{isession},-1);
%     sig_corr_all_2{isession,1} = sig_corr_all_2{isession,1}(find(sig_corr_all_2{isession,1}));
    sig_corr_all{isession,1} = tril(sig_corr_matrix{isession},-1);
    sig_corr_all{isession,1} = sig_corr_all{isession,1}(find(sig_corr_all{isession,1}));
    
    noise_corr_all_1{isession,1}  = tril(noise_corr_matrix1{isession},-1);
    noise_corr_all_1{isession,1}  = noise_corr_all_1{isession,1}(find(noise_corr_all_1{isession,1}));
    noise_corr_all_2{isession,1}  = tril(noise_corr_matrix2{isession},-1);
    noise_corr_all_2{isession,1} = noise_corr_all_2{isession,1}(find(noise_corr_all_2{isession,1})); 
    

end
Sig_corr_all = cell2mat(sig_corr_all);
Noise_corr_all_1 = cell2mat(noise_corr_all_1);
[corr_r_1,corr_p_1] = corr(Sig_corr_all,Noise_corr_all_1);
text_1 = ['r=',num2str(corr_r_1),';p=',num2str(corr_p_1)];

Sig_corr_all = cell2mat(sig_corr_all);
Noise_corr_all_2 = cell2mat(noise_corr_all_2); 
[corr_r_2,corr_p_2] = corr(Sig_corr_all,Noise_corr_all_2);
text_2 = ['r=',num2str(corr_r_2),';p=',num2str(corr_p_2)];

gray = [0.3,0.3,0.3];
figure(1)
set(figure(1),'position',[60 91 1200 500]);
popu1 = Noise_corr_all_1<0;
popu2 = Noise_corr_all_1>0;
h_subplot = tight_subplot(1,2,[0.05 0.05],[0.1 0.1],[0.08 0.03],[],[]);
set(gcf,'CurrentAxes',h_subplot(1))
plot(Sig_corr_all(popu1),Noise_corr_all_1(popu1),'o','Color',gray)
hold on
plot(Sig_corr_all(popu2),Noise_corr_all_1(popu2),'o','Color',gray)
title(text_1)
p = polyfit(Sig_corr_all, Noise_corr_all_1, 1);   % 一阶多项式
slope = p(1);
intercept = p(2);
x_fit = linspace(min(Sig_corr_all), max(Sig_corr_all), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '-', 'Color', gray, 'LineWidth', 2)
ylim([-1,1])
set(gcf,'CurrentAxes',h_subplot(2))
plot(Sig_corr_all(popu1),Noise_corr_all_2(popu1),'o','Color',gray)
hold on
plot(Sig_corr_all(popu2),Noise_corr_all_2(popu2),'o','Color',gray)
title(text_2)
p = polyfit(Sig_corr_all, Noise_corr_all_2, 1);   % 一阶多项式
slope = p(1);
intercept = p(2);
x_fit = linspace(min(Sig_corr_all), max(Sig_corr_all), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '-', 'Color', gray, 'LineWidth', 2)
ylim([-1,1])
SetFigure
figure(2)
plot(Noise_corr_all_1,Noise_corr_all_2,'o','Color',gray)
p = polyfit(Noise_corr_all_1, Noise_corr_all_2, 1);   % 一阶多项式
slope = p(1);
intercept = p(2);
x_fit = linspace(min(Noise_corr_all_1), max(Noise_corr_all_1), 100);
y_fit = polyval(p, x_fit);
hold on
plot(x_fit, y_fit, '-', 'Color', gray, 'LineWidth', 2)
xlim([-1,1])
ylim([-1,1])
[corr_r,corr_p] = corr(Noise_corr_all_1,Noise_corr_all_2);
text_3 = ['r=',num2str(corr_r),';p=',num2str(corr_p)];
title(text_3)
SetFigure
pbaspect([1 1 1])

figure(3)
% Noise_corr_all_1_sac = cell2mat(noise_corr_all_1_sac);
% Noise_corr_all_2_sac = cell2mat(noise_corr_all_2_sac);
Noise_corr_all_1_sac = vertcat(noise_corr_all_1_sac{:});
Noise_corr_all_2_sac = vertcat(noise_corr_all_2_sac{:});
plot(Noise_corr_all_1_sac,Noise_corr_all_2_sac,'o','Color',gray)
p = polyfit(Noise_corr_all_1_sac, Noise_corr_all_2_sac, 1);   % 一阶多项式
slope = p(1);
intercept = p(2);
x_fit = linspace(min(Noise_corr_all_1_sac), max(Noise_corr_all_1_sac), 100);
y_fit = polyval(p, x_fit);
hold on
plot(x_fit, y_fit, '-', 'Color', gray, 'LineWidth', 2)
xlim([-0.5,1])
ylim([-0.5,1])
[corr_r,corr_p] = corr(Noise_corr_all_1_sac,Noise_corr_all_2_sac);
text_3 = ['r=',num2str(corr_r),';p=',num2str(corr_p)];
title(text_3)
SetFigure
pbaspect([1 1 1])
