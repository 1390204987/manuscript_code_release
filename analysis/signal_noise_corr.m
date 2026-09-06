%% plot the signal corr vs noise corr of neural recorded from LA

clear; close all;

% FolderPath = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dLIP_LA\';
% FolderPath = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dVIP_LA\';
% FolderPath = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dMST_LA\';
FolderPath = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseIPS_LA\';
FolderPath_restrict = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogonIPS_LA\'; %for select neuron recorded in both task
% FolderPath = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dMST_LA\';
% FolderPath_restrict = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogon1dMST_LA\'; %for select neuron recorded in both task
% FolderPath = 'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dSTS_LA\';
% FolderPath_restrict = 'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Friday1dSTS_LA\'; %for select neuron recorded in both task
% FolderPath = 'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaymemSTS_LA\';
% FolderPath_restrict = 'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaymemSTS_LA\'; %for select neuron recorded in both task
% FolderPath = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogonmemMST_LA\';
% FolderPath_restrict = 'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogonmemMST_LA\'; %for select neuron recorded in both task
% taskname = '1Dtuning';
taskname = 'HD'; 
% taskname = 'memSac';
if strcmp(taskname,'1Dtuning')
    files = dir(fullfile([FolderPath,'*Azimuthtuning.mat']));  
elseif strcmp(taskname,'colorHD')||strcmp(taskname,'HD')
    files = dir(fullfile(FolderPath,'*PSTH.mat'));
elseif strcmp(taskname,'memSac')
    files = dir(fullfile(FolderPath,'*MemSac.mat'));
end 

taskname_restrict = 'colorHD';
% taskname_restrict = '1Dtunning';
% taskname_restrict = 'memSac';
if strcmp(taskname_restrict,'colorHD')||strcmp(taskname_restrict,'HD') 
    files_restrict = dir(fullfile(FolderPath_restrict,'*PSTH.mat'));
elseif strcmp(taskname_restrict,'1Dtunning')
    files_restrict = dir(fullfile([FolderPath_restrict,'*Azimuthtuning.mat']));  
elseif strcmp(taskname_restrict,'memSac')
    files_restrict = dir(fullfile([FolderPath_restrict,'*MemSac.mat']));  
end

for ifile = 1:length(files)
loc_ = strfind(files(ifile).name,'_');
files_name{ifile}=[files(ifile).name(1:loc_(1)-3),files(ifile).name(loc_(1):(loc_(2)-1))];
files_session{ifile} = [files(ifile).name(1:loc_(1)-3)];
end

for ifile = 1:length(files_restrict)
loc_ = strfind(files_restrict(ifile).name,'_');
files_name_restrict{ifile}=[files_restrict(ifile).name(1:loc_(1)-3),files_restrict(ifile).name(loc_(1):(loc_(2)-1))];
files_session_restrict{ifile} = [files_restrict(ifile).name(1:loc_(1)-3)];
end

[files_name_both,useless1,useless2] = intersect(files_name,files_name_restrict);
[files_session_both,useless1,useless2] = intersect(files_session,files_session_restrict);

% select_files_name = files_name;
select_files_name = files_name_both;
select_session_name = files_session_both;

unique_files_name = unique(select_files_name); % here get each unique session
unique_session_name = unique(select_session_name);
Files_session = cell(length(unique_session_name),1);

for ineuron = 1: length(unique_files_name)
neuronID = unique_files_name{ineuron};
cutline = strfind(neuronID,'_');
neuronID_session = neuronID(1:cutline-1);
neuronID_channel = neuronID(cutline+1:end);
if strcmp(taskname,'1Dtuning')
    Files = dir(fullfile([FolderPath,neuronID_session,'r','*','_',neuronID_channel,'_','*AzimuthTuning.mat']));
elseif strcmp(taskname,'colorHD')
    Files = dir(fullfile([FolderPath,neuronID_session,'r','*','_',neuronID_channel,'_','*PSTH.mat']));
elseif strcmp(taskname,'HD')
    Files = dir(fullfile([FolderPath,neuronID_session,'r','*','_',neuronID_channel,'_','*PSTH.mat']));
elseif strcmp(taskname,'memSac')
    Files = dir(fullfile([FolderPath,neuronID_session,'r','*','_',neuronID_channel,'_','*MemSac.mat']));
end

for isession = 1:length(unique_session_name)
    if strcmp(neuronID_session ,select_session_name{isession})
        Files_session{isession,1} = [Files_session{isession,1};{Files.name}];
    end
end
end

for isession = 1:length(unique_session_name)
if length(Files_session{isession,1}) < 2
    continue
end
list = [];
for ct = 1:length(Files_session{isession,1})
    list{ct,1} = load([FolderPath Files_session{isession}{ct}]); % the list contain files recorded from the same session
end

%taskname = 'colorHD';

signal_corr = signal_corr_LA(list,taskname);
noise_corr = noise_corr_LA(list,taskname);

signal_corr = tril(signal_corr,-1);
signal_corr = signal_corr(find(signal_corr));
noise_corr = tril(noise_corr,-1);
noise_corr = noise_corr(find(noise_corr));

signal_corr_all{isession,1} = signal_corr;
noise_corr_all{isession,1} = noise_corr;

% [cor_signalnoise{ineuron},Pcor_signalnoise{ineuron}] = corr(signal_corr,noise_corr);
end

Signal_corr_all = cell2mat(signal_corr_all);
Noise_corr_all = cell2mat(noise_corr_all);

figure(2)
plot(Signal_corr_all,Noise_corr_all,'b.')

[beta,p_beta] = corr(Signal_corr_all,Noise_corr_all,'rows','complete');
hline = refline(beta,0);

if  p_beta<0.05
    hline.LineWidth = 1;
end
% xlims = xlim;
% ylims = ylim;
% text(min(xlims),max(ylims)-0.1,['p=', num2str(p_beta)])
% text(min(xlims),max(ylims),['beta=', num2str(beta)])
ylim([-1,1])
xlabel(['beta=', num2str(beta),' p=', num2str(p_beta)])
title(['n=',num2str(length(Signal_corr_all))])
SetFigure
save('Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorLIP_LA\LIPcolor_signalnoisecorr.mat','Signal_corr_all','Noise_corr_all')
