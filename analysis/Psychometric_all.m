% population neuron metric
clc;clear;close all;
% FolderPath = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\Fridaycoarsecolorbehavior\';'Z:\Data\Tempo\Batch\fingerpostFriday_NP\Fridaycolorbehavior\';};
% FolderPath = {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\Drogoncoarsecolorbehavior\'};
FolderPath = {'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorNTbehavior\'};
%% load all cell color task  behavior data
if length(FolderPath)==1
    Files_psycho = dir(fullfile(FolderPath{1},'*_Psycho.mat'));
    pathfilenum = length(Files_psycho);
else
    Files_psycho = [];
    for ipath = 1:length(FolderPath)
        ipathfiles = dir(fullfile(FolderPath{ipath},'*_Psycho.mat'));
        Files_psycho = [Files_psycho;ipathfiles];
        pathfilenum(ipath) = length(ipathfiles);
    end
    pathfilenum = cumsum(pathfilenum);
end
File_work_count = 0;
for ct = 1:length(Files_psycho)
    File_work_count = File_work_count+1;
    current_path_pod = find(pathfilenum<File_work_count);
    if isempty(current_path_pod)
        current_path = 1;
    else
        current_path = max(current_path_pod) + 1;
    end
    list_psycho{File_work_count} = load([FolderPath{current_path} Files_psycho(ct).name]);
    TF = isstrprop(list_psycho{1,ct}.result.FILE,'digit');
    charloc = find(TF==0);
    %         monkeyID = str2num(list{1,ct}.result.FILE(charloc(1)+1:charloc(2)-1));
    blockID(ct) = str2double(list_psycho{1,ct}.result.FILE(charloc(2)+1:charloc(3)-1))*1000+str2double(list_psycho{1,ct}.result.FILE(charloc(3)+1:end))*10;
end

%%
blocks = length(Files_psycho)
for ct = 1:length(Files_psycho)
    heading_per_block{ct,1} = list_psycho{1,ct}.result.heading';
    motion_coherence_per_block{ct,1} = list_psycho{1,ct}.result.motion_coherence';
    motion_cal_coherence_per_block{ct,1} = list_psycho{1,ct}.result.motion_cal_coherence';
    
    outcome_per_block{ct,1} = list_psycho{1,ct}.result.outcome';
    LEFT_RF_color_block{ct,1} = list_psycho{1,ct}.result.LEFT_RF_color';
    choice_per_block{ct,1} = list_psycho{1,ct}.result.choice_per_trial';
    stim_type_block{ct,1} = list_psycho{1,ct}.result.stim_type';
    
    correct_proportion_per_block{ct,1} = list_psycho{1,ct}.result.correct_proportion;
    thresh_psy_block{ct,1} = list_psycho{1,ct}.result.Thresh_psy;
    bias_psy_block{ct,1} = list_psycho{1,ct}.result.Bias_psy;
end
heading_all = cell2mat(heading_per_block);
LEFT_RF_color_all = cell2mat(LEFT_RF_color_block);
choice_all = cell2mat(choice_per_block);

unique_heading_all = munique(heading_all);
unique_motion_coherence_all = 60;
unique_LEFT_RF_color_all = munique(LEFT_RF_color_all);

%determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)    
LEFT = 1;
RIGHT = 2;

correct_proportion = [];

if length(unique_LEFT_RF_color_all) > 1
    colors = [1,2,3];
else
    colors = [3];
end

for icolor = colors
    for icoherence = 1:length(unique_motion_coherence_all) % different coherence level
          for iheading = 1:length(unique_heading_all)
              if icolor == 3
                  trials_select =logical( (heading_all == unique_heading_all(iheading)) );
              else
                  trials_select =logical( (heading_all == unique_heading_all(iheading)) &(LEFT_RF_color_all == unique_LEFT_RF_color_all(icolor)));% & (motion_coherence==unique_motion_coherence(c))); %& data.one_targ_params==0 ) ;%changed by zn
              end
              rightward_trials = (trials_select & (choice_all == RIGHT) );
              rightward_rate = 1*sum(rightward_trials) / sum(trials_select);
              fit_data_psycho_cum{icolor,icoherence}(iheading, 1) = unique_heading_all(iheading);
              fit_data_psycho_cum{icolor,icoherence}(iheading, 2) = rightward_rate;
              fit_data_psycho_cum{icolor,icoherence}(iheading, 3) = sum(trials_select);
              %% calculate under each block for ttest
              if icolor~=3
                  for iblock = 1:length(Files_psycho)
                      trials_select_block =logical( (heading_per_block{iblock} == unique_heading_all(iheading)) ...
                          &(LEFT_RF_color_block{iblock} == unique_LEFT_RF_color_all(icolor)));%changed by zn
                      rightward_trials_block= (trials_select_block & (choice_per_block{iblock}== RIGHT) );
                      rightward_rate_block{icolor,iheading}(iblock,1) = 1*sum(rightward_trials_block) / sum(trials_select);
                  end
              end
          end
          % the correct rate does not take coherence into account,temporarily 05/29/09
          if icolor == 3
              trials_rightward = find( (heading_all > 0) & (choice_all==RIGHT) );  %&(data.one_targ_params==0)  ) ;
              trials_leftward  = find( (heading_all < 0) & (choice_all==LEFT)); %&(data.one_targ_params==0) ) ;
              trials_all = find( ((heading_all < 0)|(heading_all > 0)) ); %&(data.one_targ_params==0) ); %exclude 0 headings
          else
              trials_rightward = find( (heading_all > 0) & (choice_all==RIGHT) &(LEFT_RF_color_all == unique_LEFT_RF_color_all(icolor)));  %&(data.one_targ_params==0)  ) ;
              trials_leftward  = find( (heading_all < 0) & (choice_all==LEFT) &(LEFT_RF_color_all == unique_LEFT_RF_color_all(icolor))); %&(data.one_targ_params==0) ) ;
              trials_all = find( ((heading_all < 0)|(heading_all > 0)) &(LEFT_RF_color_all == unique_LEFT_RF_color_all(icolor))); %&(data.one_targ_params==0) ); %exclude 0 headings
          end
          correct_proportion(icolor,1) = (length(trials_rightward)+length(trials_leftward))/length(trials_all);
          
          aa = find(fit_data_psycho_cum{icolor,icoherence}(:,2)>-99); % sometime it could be NaN due to the absence of that heading conditions
          fit_valid{icolor,icoherence}(:,1) = fit_data_psycho_cum{icolor,icoherence}(aa,1);
          fit_valid{icolor,icoherence}(:,2) = fit_data_psycho_cum{icolor,icoherence}(aa,2);
          fit_valid{icolor,icoherence}(:,3) = fit_data_psycho_cum{icolor,icoherence}(aa,3);
    end
end

method = 1; % 0: Maximum likelihood; 1: Square error
tolerance = 10; 

for icolor = colors
    for icoherence = 1:length(unique_motion_coherence_all) % different coherence level      
        %   similar way to fit data
        [bb,tt] = cum_gaussfit_max1(fit_valid{icolor,icoherence},method,0);
        [bb_tol,tt_tol] = cum_gaussfit_max1(fit_valid{icolor,icoherence},method,tolerance);
        Thresh_psy{icolor,icoherence} = tt;        Bias_psy{icolor,icoherence} = bb;
        Thresh_psy_tol{icolor,icoherence} = tt_tol;    Bias_psy_tol{icolor,icoherence} = bb_tol;
        psy_perf{icolor,icoherence} =[bb,tt];
    end
end
%% ttest check if different color location cause behavior difference
if length(unique_LEFT_RF_color_all) > 1
    for iheading = 1:length(unique_heading_all)
        [~,p_heading(iheading)] = ttest(rightward_rate_block{1,iheading}(:),rightward_rate_block{2,iheading}(:));
    end
end
%%
% % plot psychometric
% % run the slide threshold over time, see whether performance fluctuate across time
% not work for coherence, temporarily 05/29/09

% plot psychometric function here
 symbo{1,1} = 'o';     fitline{1,1} = '-'; 
 symbo{2,1} = 'o';     fitline{2,1} = '-'; 
 symbo{3,1} = 'o';     fitline{3,1} = '-'; 
 color{1,1} = [205 0 0]/255;     color{1,1} = [205 0 0]/255; 
 color{2,1} = [0 139 0]/255;     color{2,1} = [0,139 0]/255;  
 color{3,1} = 'k';     color{3,1} = 'k';  
figure(2);
set(2,'Position', [200,50,700,900], 'Name', 'Heading Discrimination-Visual');
a1 = axes('position',[0.2,0.47, 0.5,0.3] );
a2 = axes('position',[0.2,0.1, 0.5,0.3] );
% fit data with cumulative gaussian and plot both raw data and fitted curve
legend_txt = [];
xi = min(unique_heading_all) : 0.1 : max(unique_heading_all);

for icolor = colors
    for icoherence = 1:length(unique_motion_coherence_all) % different coherence level
        if icolor == 3
            axes(a2);
            %                 axes(a1);
        else
            axes(a1);
        end
        plot(unique_heading_all, fit_valid{icolor,icoherence}(:,2), symbo{icolor,icoherence},  xi, cum_gaussfit(psy_perf{icolor,icoherence}, xi),  fitline{icolor,icoherence}...
            ,'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',color{icolor,icoherence},'MarkerEdgeColor',color{icolor,icoherence},'Color',color{icolor,icoherence});
        hold on
        %             if (Protocol == COARSE_HEADING_DISCRIM_COLOUR)
        xlabel('coherence');
        %             else
        %             xlabel('Heading Angles');
        %             end
        ylim([0,1]);
        ylabel('Red Choices');
        
%         set(gca, 'YTickMode','auto');
        %             set(gca, 'xTickMode','auto');
        yticks([0,0.5,1])
%         xticks([-2,-1,0,1,2])
        %             xticklabels({'-easy','-mid','-hard','0','hard','mid','easy'})
%         xticklabels({'-easy','-hard','0','hard','easy'})
        hold on;
%         legend_txt{istimtype*2-1} = [num2str(unique_stim_type_all(istimtype))];
%         legend_txt{istimtype*2} = [''];        
    end
end

% output some text of basic parameters in the figure
axes('position',[0.2,0.83, 0.6,0.15] );
xlim( [0,50] );
ylim( [2,10] );
SetFigure
text(10,8, 'u             sigma          correct rate         diffexist');
axis off;
for icolor = colors
    if icolor == 3
        for icoherence = 1:length(unique_motion_coherence_all) % different coherence level           
            text(10,7-icolor-(icoherence-1)*3,num2str(Bias_psy{icolor,icoherence}),'color',color{icolor,icoherence});
            text(20,7-icolor-(icoherence-1)*3,num2str(Thresh_psy{icolor,icoherence}),'color',color{icolor,icoherence});
            text(30,7-icolor-(icoherence-1)*3,num2str(correct_proportion(icolor)),'color',color{icolor,icoherence});
            if length(unique_LEFT_RF_color_all) > 1
                diff_exist = ~isempty(find(p_heading<0.05));
                text(40,7-icolor-(icoherence-1)*3,num2str(diff_exist),'color',color{icolor,icoherence});
            end            
        end
    else
        for icoherence = 1:length(unique_motion_coherence_all) % different coherence level            
            text(0,7-icolor-(icoherence-1)*3, num2str(unique_LEFT_RF_color_all(icolor)));  % non-microstim
            text(10,7-icolor-(icoherence-1)*3,num2str(Bias_psy{icolor,icoherence}),'color',color{icolor,icoherence});
            text(20,7-icolor-(icoherence-1)*3,num2str(Thresh_psy{icolor,icoherence}),'color',color{icolor,icoherence});
            text(30,7-icolor-(icoherence-1)*3,num2str(correct_proportion(icolor)),'color',color{icolor,icoherence});       
        end
    end
end













