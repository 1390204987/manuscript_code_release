% here plot cp value turn over under T1 on left and T1 on right condition
% temporally  by zn 2022/12/1
function [true_diverge,pvalue_CP_diverge,t_center_step] = CPturn_quantify(varargin)
if ~isempty(varargin)
    list = varargin{1};
    num_shuffles = 1000;
else
    %%
    num_shuffles = 1000;
    FolderPath = 'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Fridaycolor113up_LA\';
    neuronID = 'm24c39r3_68';
    File = dir(fullfile([FolderPath,neuronID,'*PSTH.mat']));
    % File = dir(fullfile(FolderPath,'*_PSTH.mat'));
    list = load([FolderPath File.name]);
end
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
% spike_hist = list.result.spike_hist{1,1};
spike_aligned = list.result.spike_aligned;
align_offsets_others = list.result.align_offsets_others;
spike_rates = list.result.spike_rates;
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
time_window2(1,1) = find(spike_aligned{2,1}(:) > stimulus_on,1);
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

%%
% zscore
% shift window, currently using 100ms shift with 1000ms width
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
        z_heading_dist_temp = z_heading_dist - mean(z_heading_dist,1)/std(z_heading_dist(:));

        Zspike_heading_pertrial(select,:) = z_heading_dist_temp;
        z_heading_dist_temp = [];
    end   
end
%%

num_trials = length(choice_per_trial); % 获取 choice_per_trial 的长度


% 使用 arrayfun 生成 1000 次随机打乱
shuffled_indices = arrayfun(@(x) randperm(num_trials), 1:num_shuffles, 'UniformOutput', false);

% 如果需要打乱后的 choice_per_trial 数据
shuffled_T1loc = cellfun(@(idx) T1loc_per_trial(idx), shuffled_indices, 'UniformOutput', false);
T1loc = {T1loc_per_trial,shuffled_T1loc{:}}; % the first column is true lable the others is shuffled label
Zresp_left_choice_T1L = cellfun(@(x) Zspike_heading_pertrial((choice_per_trial==-1)&(sign(x)==-1),:),T1loc, 'UniformOutput', false);
Zresp_right_choice_T1L = cellfun(@(x) Zspike_heading_pertrial((choice_per_trial==1)&(sign(x)==-1),:),T1loc,'UniformOutput', false);

Zresp_left_choice_T1R = cellfun(@(x) Zspike_heading_pertrial((choice_per_trial==-1)&(sign(x)==1),:),T1loc,'UniformOutput', false);
Zresp_right_choice_T1R = cellfun(@(x) Zspike_heading_pertrial((choice_per_trial==1)&(sign(x)==1),:),T1loc,'UniformOutput', false);

%% T test evaluate choice selectivity
for s = 1:length(t_center_step{1})
%     if length(Zresp_left_choice_T1L(:,s)) >3 && length(Zresp_right_choice_T1L(:,s))>3
        choice_roc_T1L(s,:)= cellfun(@(x,y) rocN(x(:,s),y(:,s)),Zresp_left_choice_T1L,Zresp_right_choice_T1L);  %set the left as x axis
%     else
%         [choice_roc_T1L(s,1),perm] = NaN;
%     end
    
%     if length(Zresp_left_choice_T1R(:,s)) >3 && length(Zresp_right_choice_T1R(:,s))>3
        choice_roc_T1R(s,:)= cellfun(@(x,y) rocN(x(:,s),y(:,s)),Zresp_left_choice_T1R,Zresp_right_choice_T1R);  %set the left as x axis
%     else
%         [choice_roc_T1R(s,1),perm] = NaN;
%     end
end
[rr,pp] = corrcoef(relativeheading_per_trial, spike_rates');

line(1,1) = rr(1,2);
line(1,2) = pp(1,2);

if  logical(line(1,1)>0)%&logical(~isnan(choice_roc_T1L))
    CP_T1L = 1-choice_roc_T1L;
    CP_T1R = 1-choice_roc_T1R;
else
    CP_T1L = choice_roc_T1L;
    CP_T1R = choice_roc_T1R;
end
CP_T1L = choice_roc_T1L;
CP_T1R = choice_roc_T1R;

abs_CP_diverge = abs(CP_T1R-CP_T1L);
pvalue_CP_diverge = sum(abs_CP_diverge(:,2:end)>abs_CP_diverge(:,1),2)/num_shuffles;
mean_shuffle_diverge = mean(abs_CP_diverge(:,2:end),2);
true_diverge = abs_CP_diverge(:,1);
%% figure plot
if isempty(varargin)
    diverge_color = [0,0,0];
    shuffle_color = [0.5,0.5,0.5];
    figure(1); clf;
    set(figure(1),'position',[60 91 600 200])
    hold on
    h_legend(1) = plot(t_center_step{1},abs_CP_diverge(:,1),'linewidth',1,'Color',diverge_color,'marker','o','markersize',10,'MarkerEdgeColor',diverge_color,'linestyle','-');
    plot(t_center_step{1}(CP_T1L(:,2)<0.05),CP_T1L(CP_T1L(:,2)<0.05,1),...
        ['o'], 'markerfacecolor',diverge_color,'markersize',10,'MarkerEdgeColor',diverge_color);   % Sig. CP
    
    h_legend(2) = plot(t_center_step{1},mean_shuffle_diverge,'linewidth',1,'Color',shuffle_color,'marker','o','markersize',10,'MarkerEdgeColor',shuffle_color,'linestyle','-');
    plot(t_center_step{1}(pvalue_CP_diverge<0.05),abs_CP_diverge(pvalue_CP_diverge<0.05,1),...
        ['o'], 'markerfacecolor',shuffle_color,'markersize',10,'MarkerEdgeColor',diverge_color);   % Sig. CP
    %  plot([t_center_step{1}(1),t_center_step(end)],[0.5,0.5],'--','linewidth',0.8,'Color','k')
    xlim([0,1500])
     
    %  ylim([0.15, 0.85])
    %  title(num2str(cellID))
    SetFigure
end
end

