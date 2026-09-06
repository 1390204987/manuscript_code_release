% here plot cp value turn over under T1 on left and T1 on right condition
% temporally  by zn 2022/12/1
close all;
clear; clc
%%
FolderPath = 'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTIPS_NP\';
% FolderPath = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseIPS_LA\';
neuronID = 'm24c115r3_293';
% neuronID = 'm24c115r5_293';
File = dir(fullfile([FolderPath,neuronID,'*PSTH.mat']));
% File = dir(fullfile(FolderPath,'*_PSTH.mat'));
for ifile = 1:length(File)
list = load([FolderPath File(ifile).name]);

cell_channel = list.result.SpikeChan;
TF = isstrprop(list.result.FILE,'digit');
charloc = find(TF==0);
monkeyID = str2num(list.result.FILE(charloc(1)+1:charloc(2)-1));
cellnum = str2double(list.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;
blocknum = str2double(list.result.FILE(charloc(3)+1:end));
cellID = str2double(list.result.FILE(charloc(2)+1:charloc(3)-1))*10^5+cell_channel;   %+blocknum(ct)/10;

heading_per_trial = (list.result.heading_per_trial)';
coherence_per_trial = (list.result.coherence_per_trial)';
T1loc_per_trial = (list.result.targetloc_per_trial)';
choice_per_trial = (list.result.choice_per_trial)';
sac_per_trial = sign(T1loc_per_trial).*choice_per_trial; %
% spike_hist = list.result.spike_hist{1,1};
spike_aligned = list.result.spike_aligned;
align_offsets_others = list.result.align_offsets_others;
% spike_rates = list.result.spike_rates;
ts_per_trial = list.result.t_centers{1,1}';
align_markers = (list.result.align_markers);
%time window
stim_dur = align_offsets_others{1,2}(1,2) - align_offsets_others{1,2}(1,1);
stimulus_on = 0;
stimulus_end = 0+stim_dur;
%method1 calcalate with psth data
time_window1(1,1) = find(ts_per_trial(:) > stimulus_on-800,1);
time_window1(1,2) = find(ts_per_trial(:) < stimulus_end,1,'last');
time_window1(1,3) = size(ts_per_trial,1);
%method2 calculate with 0/1 spike data
time_window2(1,1) = find(spike_aligned{2,1}(:) > stimulus_on+500,1);
time_window2(1,2) = find(spike_aligned{2,1}(:) < stimulus_end,1,'last');
time_window2(1,3) = size(spike_aligned{2,1},2);

if length(unique(coherence_per_trial))==1
    relativeheading_per_trial = heading_per_trial;
else
    relativeheading_per_trial = sign(heading_per_trial).*coherence_per_trial;
end
unirelativeheading = unique(relativeheading_per_trial);
unichoice = unique(choice_per_trial);
% spikehist_per_trial = spike_hist;
spike_rates = sum(spike_aligned{1,1}(:,time_window2(1,1):time_window2(1,2)),2)/((time_window2(1,2)-time_window2(1,1))/1000);

%%
% zscore
% % shift window, currently using 100ms shift with 1000ms width
% window_interval = 10;
% step = 10;
% endloop = floor((time_window1(1,3)-time_window1(1,1)-window_interval/2)/step);
% t_center_step = ts_per_trial(time_window1(1,1):step:time_window1(1,1)+step*endloop);
% for iheading = 1:length(unirelativeheading)
%     select = (relativeheading_per_trial == unirelativeheading(iheading));
%     z_heading_dist = spikehist_per_trial(select,:);
%     for s = 1:endloop+1
%         vertical = time_window1(1,1)-0.5*window_interval+step*(s-1): time_window1(1,1)+ 0.5*window_interval+step*(s-1) ;
%         z_heading_dist_step(:,s) = (mean(z_heading_dist(:,vertical),2) - mean(mean(z_heading_dist(:,vertical),2)))/std(mean(z_heading_dist(:,vertical),2),1);
%     end
%     Zspike_heading_pertrial(select,:) = z_heading_dist_step;
%     z_heading_dist_step = [];
% end

% calculate spike_hist with specific time bin
% Anne Churchland NN, bin = 10 ms with Gaussian filter (sigma = 50 ms)
binSize_rate = 100;  % in ms
stepSize_rate = 100; % in ms
smoothFactor = 125; % in ms !!
spike_timeWin = 1;
for j = 1   % For each desired marker
    t_center_step{j} = align_markers{j,2} + binSize_rate/2 : stepSize_rate : align_markers{j,3} - binSize_rate/2; % Centers of PSTH time windows
    
    spike_hist{j} = zeros(size(spike_aligned{1,j},1),length(t_center_step{j})); % Preallocation
    for k = 1:length(t_center_step{j})
        winBeg = ceil(((k-1) * stepSize_rate) / spike_timeWin) + 1;
        winEnd = ceil(((k-1) * stepSize_rate + binSize_rate) / spike_timeWin) + 1;
        spike_hist{j}(:,k) = sum(spike_aligned{1,j}(: , winBeg:winEnd),2) / binSize_rate*1000 ;  % in Hz
    end
    
    % --- Smoothing ---
    % Anne Churchland NN, bin = 10 ms with Gaussian filter (sigma = 50 ms)
    % Note that this only influence PSTH calculation (spike_hist), not CP
    if smoothFactor > 0
        for i = 1:size(spike_hist{j},1) % Each trial
            spike_hist{j}(i,:) =  GaussSmooth(t_center_step{j},spike_hist{j}(i,:),smoothFactor);
        end
    end
    for iheading = 1:length(unirelativeheading)
        select = (relativeheading_per_trial == unirelativeheading(iheading));
        z_heading_dist =  spike_hist{j}(select,:);
        z_heading_dist_temp = (z_heading_dist - mean(z_heading_dist,1))/std(z_heading_dist(:));

        Zspike_heading_pertrial(select,:) = z_heading_dist_temp;
%         Zspike_heading_pertrial(select,:) = z_heading_dist;
        z_heading_dist_temp = [];
    end   
end
%%

Zresp_left_choice_T1L = Zspike_heading_pertrial((choice_per_trial==-1)&(sign(T1loc_per_trial)==-1),:);
Zresp_right_choice_T1L = Zspike_heading_pertrial((choice_per_trial==1)&(sign(T1loc_per_trial)==-1),:);

Zresp_left_choice_T1R = Zspike_heading_pertrial((choice_per_trial==-1)&(sign(T1loc_per_trial)==1),:);
Zresp_right_choice_T1R = Zspike_heading_pertrial((choice_per_trial==1)&(sign(T1loc_per_trial)==1),:);

%% T test evaluate choice selectivity
for s = 1:length(t_center_step{1})
    if length(Zresp_left_choice_T1L(:,s)) >3 && length(Zresp_right_choice_T1L(:,s))>3
        [choice_roc_T1L(s,1),~,perm] = rocN(Zresp_left_choice_T1L(:,s),Zresp_right_choice_T1L(:,s),100,1000);  %set the left as x axis
        choice_roc_T1L(s,2) = perm.pValue;
%         [~,choice_roc_T1L(s,2)] = ttest2(Zresp_left_choice_T1L(:,s),Zresp_right_choice_T1L(:,s)); % Use ttest2 instead of permutation for speed (the outcomes are quite similar). HH20140601
    else
        [choice_roc_T1L(s,1),perm] = NaN;
    end
    
    if length(Zresp_left_choice_T1R(:,s)) >3 && length(Zresp_right_choice_T1R(:,s))>3
        [choice_roc_T1R(s,1),~,perm] = rocN(Zresp_left_choice_T1R(:,s),Zresp_right_choice_T1R(:,s),100,1000);  %set the left as x axis
        choice_roc_T1R(s,2) = perm.pValue;
%         [~,choice_roc_T1R(s,2)] = ttest2(Zresp_left_choice_T1R(:,s),Zresp_right_choice_T1R(:,s)); % Use ttest2 instead of permutation for speed (the outcomes are quite similar). HH20140601
    else
        [choice_roc_T1R(s,1),perm] = NaN;
    end    
end
[rr,pp] = corrcoef(relativeheading_per_trial, spike_rates');

line(1,1) = rr(1,2);
line(1,2) = pp(1,2);

if  logical(line(1,1)>0)&logical(~isnan(choice_roc_T1L))
    CP_T1L(:,1) = 1-choice_roc_T1L(:,1);
    CP_T1R(:,1) = 1-choice_roc_T1R(:,1);
else
    CP_T1L(:,1) = choice_roc_T1L(:,1);
    CP_T1R(:,1) = choice_roc_T1R(:,1);
end
CP_T1L(:,2) = choice_roc_T1L(:,2);
CP_T1R(:,2) = choice_roc_T1R(:,2);
%% figure plot
T1R_color = [0 139 0]/255;
T1L_color = [205 0 0]/255;   
figure(1); clf;
set(figure(1),'position',[60 91 600 300])
hold on
xrange = t_center_step{1}<1600;
h_legend(1) = plot(t_center_step{1}(xrange),CP_T1L(xrange,1),'linewidth',1,'Color',T1L_color,'marker','o','markersize',10,'MarkerEdgeColor',T1L_color,'linestyle','-');
plot(t_center_step{1}(CP_T1L(xrange,2)<0.05),CP_T1L(CP_T1L(xrange,2)<0.05,1),...
      ['o'], 'markerfacecolor',T1L_color,'markersize',10,'MarkerEdgeColor',T1L_color);   % Sig. CP
  
h_legend(2) = plot(t_center_step{1}(xrange),CP_T1R(xrange,1),'linewidth',1,'Color',T1R_color,'marker','o','markersize',10,'MarkerEdgeColor',T1R_color,'linestyle','-');
plot(t_center_step{1}(CP_T1R(xrange,2)<0.05),CP_T1R(CP_T1R(xrange,2)<0.05,1),...
      ['o'], 'markerfacecolor',T1R_color,'markersize',10,'MarkerEdgeColor',T1R_color);   % Sig. CP
 plot(t_center_step{1}(xrange),0.5*ones(size(t_center_step{1}(xrange))),'--','linewidth',0.5,'Color','k')
 
 % get p value
 [~,pvalue_CP_diverge,~] = CPturn_quantify(list);
 sigdiff_loc = find(pvalue_CP_diverge(xrange)<0.05);
%  plot(t_center_step{1}(sigdiff_loc),ones(size(sigdiff_loc))*0.5,'s','Color',[0,0,1],'MarkerFaceColor',[0,0,1],'MarkerSize',10);
plot_sigbar(t_center_step{1}(sigdiff_loc),ones(size(sigdiff_loc))*0.8,t_center_step{1}(2)-t_center_step{1}(1),[150 150 150]/255)
%  ylim([0.15, 0.85])
 ylim([0 1])
 xlim([0, 1500])
 title(num2str(cellID))
 xticks([0, 500,1000, 1500])
 yticks([0.1,0.5,0.9])
SetFigure
% figure_name = [num2str(cellID),'.png'];
% cd Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridayIPSCPturn\
% cd Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogonIPSCPturn\
% saveas(gcf,figure_name)
% legend_text = {'choose rightward heading saccade leftword','choose rightward heading saccade rightward'};
% legend(h_legend,legend_text)
end