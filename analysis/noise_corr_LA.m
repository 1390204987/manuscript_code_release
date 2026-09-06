%% 2022/9/18 by zn for     noise correlation analysis of neuroal recorded by linear array    
% calculate noise correlation between two units
function [corr_matrix] = noise_corr_LA(varargin)
% varargin: list of neural data

if isempty(varargin)
%     FolderPath = 'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dIPS_LA\';
    FolderPath = 'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaymemIPS_LA\';
    neuronID = 'm24c39r2';
%     neuronID = 'm24c39r1';
%     taskname = '1Dtuning';
    %taskname = 'colorHD';
    taskname = 'memSac';
    if strcmp(taskname,'1Dtuning')
        Files = dir(fullfile([FolderPath,neuronID,'*AzimuthTuning.mat']));
    elseif strcmp(taskname,'colorHD')
        Files = dir(fullfile([FolderPath,neuronID,'*PSTH.mat']));
    elseif strcmp(taskname,'memSac')
        Files = dir(fullfile([FolderPath,neuronID,'*MemSac.mat']));
    end
    for ct = 1:length(Files)
        list{ct,1} = load([FolderPath Files(ct).name]);
    end
    
%     list{1,1} = load('Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorMST_LA\m17c228r3_5_coarsecolor_LA_PSTH.mat');
%     list{2,1} = load('Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorMST_LA\m17c228r3_6_coarsecolor_LA_PSTH.mat');
else
    list = varargin{1,1};
    taskname = varargin{1,2};
end


spike_per_cell = cell(length(list),1);
for f = 1:length(list)
    if strcmp(taskname,'1Dtuning')
        spike_per_cell{f,1} = list{f,1}.result.all_spike_rates;
%         azimuth_per_trial{f,1} = list{f,1}.result.all_azimuth;
        % ------------- zscore for per cell----------------
        % zscore for azimuth
%         unique_azimuth = unique(azimuth_per_trial{f,1});
%         for iazimuth = 1:length(unique_azimuth)
%             select = azimuth_per_trial{f,1}==unique_azimuth(iazimuth);
%             if std(spike_per_cell{f,1}(select))==0
%                 zheading = (spike_per_cell{f,1}(select)-mean(spike_per_cell{f,1}(select)));
%             else
%                 zheading = (spike_per_cell{f,1}(select)-mean(spike_per_cell{f,1}(select)))/std(spike_per_cell{f,1}(select));
%             end
%             zheading_spike_per_cell{f,1}(select) = zheading;
%         end

    elseif strcmp(taskname,'colorHD')||strcmp(taskname,'HD')
    
        spike_hist{f,1} = [{list{f,1}.result.spike_hist{1,1}},list{f,1}.result.t_centers{1,1}'];
        spike_aligned{f,1} = [{list{f,1}.result.spike_aligned{1,1}},{list{f,1}.result.spike_aligned{2,1}}];
        
%         heading_per_trial{f,1} = (list{f,1}.result.heading_per_trial)';
%         coherence_per_trial{f,1} = (list{f,1}.result.coherence_per_trial)';
%         if length(unique(coherence_per_trial{f,1})) >1
%             relativeheading_per_trial{f,1} = sign(heading_per_trial{f,1}).*coherence_per_trial{f,1};
%         else
%             relativeheading_per_trial{f,1} = heading_per_trial{f,1};
%         end
%         unique_relativeheading{f,1} = unique(relativeheading_per_trial{f,1});
%         choice_per_trial{f,1} = (list{f,1}.result.choice_per_trial)';
%         T1loc_per_trial{f,1} = (list{f,1}.result.targetloc_per_trial)';
%         saccade_per_trial{f,1} = ((T1loc_per_trial{f,1})'.*list{f,1}.result.choice_per_trial)';
%         saccade_per_trial{f,1} = sign(saccade_per_trial{f,1});
% 
%         % check whether the T1 loc change will affect signal noise corr
% %         select_condition = sign(T1loc_per_trial{f,1})==-1; % select T1 on left
%         select_condition = sign(T1loc_per_trial{f,1})== 1; % select T1 on right
        
        %for time window selection
        stim_dur(f,1) = list{f,1}.result.align_offsets_others{1,2}(1,2) - list{1,1}.result.align_offsets_others{1,2}(1,1);
        stimulus_on(f,1) = 0;
        stimulus_end(f,1) = 0+stim_dur(f,1);
        
        %method1 calcalate with psth data
        hist_time_window1(f,1) = find(spike_hist{f,1}{1,2}(:) > stimulus_on(f,1)+500,1);
        hist_time_window1(f,2) = find(spike_hist{f,1}{1,2}(:) < stimulus_end(f,1)+1500,1,'last');
        %method2 calculate with 0/1 spike data
        bin_time_window1(f,1) = find(spike_aligned{f,1}{1,2}(:) > stimulus_on(f,1)+0,1);
        bin_time_window1(f,2) = find(spike_aligned{f,1}{1,2}(:) < stimulus_end(f,1)+1500,1,'last');

        spike_per_cell{f,1} = mean(spike_hist{f,1}{1,1}(:,hist_time_window1(f,1):hist_time_window1(f,2)),2)';
        % ------------- zscore for per cell----------------
        % zscore for heading
%         unique_heading = unique(relativeheading_per_trial{f,1});
%         for iheading = 1:length(unique_heading)
%             select = (relativeheading_per_trial{f,1}==unique_heading(iheading))&select_condition;
%             if std(spike_per_cell{f,1}(select))==0
%                 zheading = (spike_per_cell{f,1}(select)-mean(spike_per_cell{f,1}(select)));
%             else
%                 zheading = (spike_per_cell{f,1}(select)-mean(spike_per_cell{f,1}(select)))/std(spike_per_cell{f,1}(select));
%             end
%             zheading_spike_per_cell{f,1}(select) = zheading;
%         end
    elseif strcmp(taskname,'memSac')
        
        spike_per_cell{f,1} = reshape(list{f,1}.result.resp_trial{3},1,[]);
        
        
    end

end
spike_percell_mat = cell2mat(spike_per_cell);
zspike_per_cell = spike_percell_mat-mean(spike_percell_mat);
% figure(1)
% var1 = zheading_spike_per_cell{1,1};
% var2 = zheading_spike_per_cell{2,1};
% plot(var1,var2,'.');
% 
% [noisecorr ,p_noisecorr] = corr(var1,var2,'type','Pearson');
% intercept = 0.5-0.5*noisecorr;
% hline = refline(noisecorr,intercept);
% text(0.6,0.6,['p=', num2str(p_noisecorr)])


% eachlength = cellfun(@length,zspike_per_cell);
% units_fir_mat = NaN(max(eachlength),length(eachlength));
% for ineuron = 1:length(zspike_per_cell);
%     trialsnum = length(zspike_per_cell{ineuron,1});
%     units_fir_mat(1:trialsnum,ineuron) = zspike_per_cell{ineuron,1}(:);
% end
% units_fir_mat = (cell2mat(zheading_spike_per_cell))';
corr_matrix = corr(zspike_per_cell',zspike_per_cell','Rows','complete');
end