%%  Seq PCA
function Seq_PCA( TDR_name, data_name,plot_fig,save_path,save_name)
if nargin == 0
    TDR_name = 'D_IPScolor_all_mTDR.mat';
    data_name = 'D_IPScolor_all.mat';
end
Bhat=[];Bhat0=[];lambhat=[];r=[];Shat=[];What=[];
cellID_pick=[];hk=[];t_centers=[];Xmat=[];Ymat=[];
load(TDR_name)
load(data_name)
if plot_fig
    figurepath = './D_IPScolor/';
end
% filename = 'F_colorvariable_shuffleall_STS.mat';
time = t_centers{:};
P = 7;
% Remove condition-independent component from design matrix for MML
% estimation since we assume full rank for this component for all models
Xb = Xmat;     Xb(:,P) = [];
[neuro_num,T,~] = size(Ymat);
 for p = 1:4
    startlen = 1;
    startpoint = 1;
    Xb_reshape = reshape(Xb(:,p),1,1,[]);% transfer b to (1,1,Nmax)
    effectW = What{p}(:,:);
    effectS = Shat{p}(:,:);
    Bptrue{p} = effectW*effectS;
    [psth,psth_correct,psth_error,unique_cons,unique_cohe] = get_PSTH(Xmat,Ymat,p,hk);
    plotB = Bptrue{p}; % Bhat & Bptrue not perfect match
    % plotB = Bhat{p}; % Bhat & Bptrue not perfect match
    
    color_{1,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
    color_{2,1} = [0,191,255;155,48,255]/255; %color
    color_{3,1} = [139,10,80;72,61,139]/255; %choice
    color_{4,1} = [0,0,139;205,133,0]/255; %saccade
    color_{5,1} = [0,0,139;205,133,0]/255; %last saccade
    color_{6,1} = [0.1,0.1,0.1;]; %coherence
    color_{7,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
    linesize = 2;
    color = color_{p,1};
    variable_name={'heading';'color';'choice';'saccade';'last_sac';'coherence'};
    topEV{1,1} = gettopEV(Xb_reshape,plotB,startlen,startpoint);   
    stupid_time = find(time>-1000&time<50);
    startseqtime = time(find(topEV{1,1}==min(topEV{1,1}(stupid_time))));
    [seqPC(:,1),seqtime1]=getseqPC1(topEV{1,1},time,startseqtime,Xb_reshape,plotB);
    seqtime(1,1) = seqtime1;
    if ~isnan(seqtime1)
        subperiod1_B = plotB-seqPC(:,1)*seqPC(:,1)'*plotB;
        topEV{2,1} = gettopEV(Xb_reshape,subperiod1_B,startlen,startpoint);
        startseqtime = seqtime1;
        [seqPC(:,2),seqtime2]=getseqPC1(topEV{2,1},time,startseqtime,Xb_reshape,subperiod1_B);
    else
        seqtime2 = nan;
        seqPC(:,2) = nan(neuro_num,1);
    end
    seqtime(2,1) = seqtime2;
    if ~isnan(seqtime2)
        subperiod2_B = subperiod1_B-seqPC(:,2)*seqPC(:,2)'*subperiod1_B;
        topEV{3,1} = gettopEV(Xb_reshape,subperiod2_B,startlen,startpoint);
        startseqtime = seqtime2;
        [seqPC(:,3),seqtime3]=getseqPC1(topEV{3,1},time,startseqtime,Xb_reshape,subperiod2_B);
    else
        seqtime3 = nan;       
        seqPC(:,3) = nan(neuro_num,1);
        seqPC(:,2) = nan(neuro_num,1);
    end
    if isnan(seqtime3)
        seqPC(:,3) = nan(neuro_num,1);
    end
    seqtime(3,1) = seqtime3;
    if plot_fig
        if ~isnan(seqtime1)
            %plot psth projection on seqPC(:,1)
            %     [psth_heading,unique_cons,unique_cohe] = get_PSTH(Xmat,Ymat,1,hk);
            set(figure,'position',[60 91 900 300]);
            plotpsth1D(time,psth,unique_cons,unique_cohe,seqPC(:,1),color,linesize)
            ylims = ylim;
            plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
            plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
            plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
            ax = gca;
            ylabel([variable_name{p},'axis1'])
            SetFigure
            figname = [variable_name{p},'_seq1psth.png'];  % 文件名和扩展名
            saveas(gcf, fullfile(figurepath, figname));
            if p~=6 % coherence variable
                %plot divergence projection on seqPC(:,1)
                set(figure,'position',[60 91 900 300]);
    %             seqPC1_diverge(:,:,p) = plot_divergence(time,psth,unique_cons,unique_cohe,seqPC(:,1),color,linesize);
                plotDV(time,psth,unique_cons,unique_cohe,seqPC(:,1),color,linesize);
                figname = [variable_name{p},'_seq1DV.png'];
                ylims = ylim;
                plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
                plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
                plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
                ax = gca;
                ylabel([variable_name{p},'axis1'])
                SetFigure
                saveas(gcf, fullfile(figurepath, figname));
                if p==1
                    %specific for plot heading divergence
                    set(figure,'position',[60 91 900 300]);
                    [psth_2v,heading_cons,choice_cons,unique_cohe]  = get_2v_PSTH(Xmat,Ymat,hk);
                    psth_2v = permute(psth_2v,[1,2,4,3]);
                    plotdivergence(time,psth_2v,heading_cons,choice_cons,unique_cohe,seqPC(:,1),color,linesize);
                    figname = [variable_name{p},'_seq1divergence.png'];
                    ylims = ylim;
                    plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
                    plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
                    plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
                    ax = gca;
                    ylabel([variable_name{p},'axis1'])
                    SetFigure
                    saveas(gcf, fullfile(figurepath, figname));
                elseif p==3
                    %specific for plot choice divergence
                    set(figure,'position',[60 91 900 300]);
                    [psth_2v,heading_cons,choice_cons,unique_cohe]  = get_2v_PSTH(Xmat,Ymat,hk);
                    plotdivergence(time,psth_2v,choice_cons,heading_cons,unique_cohe,seqPC(:,1),color,linesize);
                    figname = [variable_name{p},'_seq1divergence.png'];
                    ylims = ylim;
                    plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
                    plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
                    plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
                    ax = gca;
                    ylabel([variable_name{p},'axis1'])
                    SetFigure
                    saveas(gcf, fullfile(figurepath, figname));
                            % specific for plot  divergence
                else
                    set(figure,'position',[60 91 900 300]);
                    control_cons = [];
                    plotdivergence(time,psth,unique_cons,control_cons,unique_cohe,seqPC(:,1),color,linesize);
                    figname = [variable_name{p},'_seq1divergence.png'];
                    ylims = ylim;
                    plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
                    plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
                    plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
                    ax = gca;
                    ylabel([variable_name{p},'axis1'])
                    SetFigure
                    saveas(gcf, fullfile(figurepath, figname));
                end

    %             set(figure,'position',[60 91 900 300]);
    %             plotDV(time,psth_correct,unique_cons,unique_cohe,seqPC(:,1),color,linesize);
    %             figname = [variable_name{p},'_seq1correctDV.png'];
    %             ylims = ylim;
    %             plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
    %             plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
    %             plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
    %             ax = gca;
    %             ylabel([variable_name{p},'axis1'])
    %             SetFigure
    %             saveas(gcf, fullfile(figurepath, figname));
    %             
    %             set(figure,'position',[60 91 900 300]);
    %             plotDV(time,psth_error,unique_cons,unique_cohe,seqPC(:,1),color,linesize);
    %             ylims = ylim;
    %             plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
    %             plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
    %             plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
    %             ax = gca;
    %             ylabel([variable_name{p},'axis1'])
    %             SetFigure
    %             figname = [variable_name{p},'_seq1errorDV.png'];
    %             saveas(gcf, fullfile(figurepath, figname));
    %         
            end
        end
        if ~isnan(seqtime2)
            %plot psth projection on seqPC(:,2)
            set(figure,'position',[60 91 900 300]);
            plotpsth1D(time,psth,unique_cons,unique_cohe,seqPC(:,2),color,linesize)
            ylims = ylim;
            plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
            plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
            plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
            ax = gca;
            ylabel([variable_name{p},'axis2'])
            SetFigure
            figname = [variable_name{p},'_seq2psth.png'];  % 文件名和扩展名
            saveas(gcf, fullfile(figurepath, figname));
            if p~=6
                %plot divergence projection on seqPC(:,2)
                set(figure,'position',[60 91 900 300]);
    %             seqPC2_diverge(:,:,p)=plot_divergence(time,psth,unique_cons,unique_cohe,seqPC(:,2),color,linesize);
                plotDV(time,psth,unique_cons,unique_cohe,seqPC(:,2),color,linesize);
                ylims = ylim;
                plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
                plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
                plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
                ax = gca;
                ylabel([variable_name{p},'axis2'])
                SetFigure
                figname = [variable_name{p},'_seq2diverge.png'];  % 文件名和扩展名
                saveas(gcf, fullfile(figurepath, figname));
            end
        end
        if ~isnan(seqtime3)
            %plot psth projection on seqPC(:,3)
            set(figure,'position',[60 91 900 300]);
            plotpsth1D(time,psth,unique_cons,unique_cohe,seqPC(:,3),color,linesize)
            plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
            plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
            plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
            ax = gca;
            ylabel([variable_name{p},'axis3'])
            SetFigure
            figname = [variable_name{p},'_seq3psth.png'];  % 文件名和扩展名
            saveas(gcf, fullfile(figurepath, figname));
            if p~=6
                %plot divergence projection on seqPC(:,3)
                set(figure,'position',[60 91 900 300]);
    %             seqPC3_diverge(:,:,p)=plot_divergence(time,psth,unique_cons,unique_cohe,seqPC(:,3),color,linesize);
                plotDV(time,psth,unique_cons,unique_cohe,seqPC(:,3),color,linesize);
                ylims = ylim;
                plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
                plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
                plot([seqtime3,seqtime3],ylims,'k','LineWidth',2)
                ax = gca;
                ylabel([variable_name{p},'axis3'])
                figname = [variable_name{p},'_seq3diverge.png'];  % 文件名和扩展名
                saveas(gcf, fullfile(figurepath, figname));
            end
        end
    end
    variablePC{p} = [seqPC(:,1),seqPC(:,2),seqPC(:,3)];
    
    % plot top EV 
    set(figure,'position',[60 91 900 300]);
    axis off
    title_infor = ['variable P  = ',num2str(p)];
    title(title_infor)
    h_subplot = tight_subplot(1,3,[0.01 0.02],[0.1 0.1],[0.08 0.03],[],[]);    
    for i = 1:3
        if i==1||(i>1&&~isnan(seqtime(i-1)))
            set(gcf,'CurrentAxes',h_subplot(i))
            hold on
            plot(time,topEV{i,1},".-")
            ylabel(['topEV',num2str(i)])
            xlabel("time")
            ylims = ylim;
            plot([seqtime1,seqtime1],ylims)
            plot([seqtime2,seqtime2],ylims)
            plot([seqtime3,seqtime3],ylims)
        end
    end
    SetFigure
    figname = [variable_name{p},'_topEV.png'];  % 文件名和扩展名
%     saveas(gcf, fullfile(figurepath, figname));
    if plot_fig
        % plot B projection
        set(figure,'position',[60 91 900 300]);
        axis off
        title_infor = ['variable P  = ',num2str(p)];
        title(title_infor)
        h_subplot = tight_subplot(1,3,[0.01 0.02],[0.1 0.1],[0.08 0.03],[],[]);    
        for i = 1:3
            if i==1||(i>1&&~isnan(seqtime(i-1)))
                set(gcf,'CurrentAxes',h_subplot(i))
                hold on
                B_project = seqPC(:,i)'*plotB;
                plot(time,B_project)
                ylims = ylim;
                plot([seqtime1,seqtime1],ylims)
                plot([seqtime2,seqtime2],ylims)
                plot([seqtime3,seqtime3],ylims)
            end
        end
        SetFigure
        figname = [variable_name{p},'_Bprojection.png'];  % 文件名和扩展名
        saveas(gcf, fullfile(figurepath, figname));

        figure
        if ~isnan(seqtime1)&&~isnan(seqtime2)
            plot2D(time,psth,unique_cons,unique_cohe,seqPC(:,1),seqPC(:,2),color,linesize)
            ax = gca;
            ax.YTick = [];
            yticklabels([]);
            ax.XTick = [];
            xticklabels([]);
            xlabel([variable_name{p},'axis1'])
            ylabel([variable_name{p},'axis2'])
            SetFigure
            figname = [variable_name{p},'2Dpsth.png'];  % 文件名和扩展名
            saveas(gcf, fullfile(figurepath, figname));
        end
    end
%     figure
%     plot3D(time,psth,unique_cons,unique_cohe,seqPC(:,1),seqPC(:,2),seqPC(:,3),color,linesize)    
%     SetFigure
 end
%plot each variable divergence together
%  figure
%  plotdivergeall(time,seqPC1_diverge) 
%  SetFigure
 save([save_path,save_name],'variablePC','cellID_pick')

%%
% load('D_variablePC_ips.mat')
% close all
% index = randperm(length(variablePC{1}(:,1))); 
% V1 = variablePC{1}(:,1);
% V2 = variablePC{3}(:,1);
% [~,neural_id] = sort(V1);
% set(figure,'position',[60 91 600 100]);
% SetFigure
% plot_weight(V1,V2,neural_id);
% 
% plot_weight_corr(V1,V2,1,'o','MarkerEdgeColor','k','MarkerFaceColor','k')

function plot_weight(V1,V2,neural_id)
blue = [0 0 255]/255;
green = [0 255 0]/255;
bar(1:length(neural_id),V1(neural_id,1),'BarWidth',1,'FaceColor',[255/255 165/255 0/255],'EdgeColor',[255/255 165/255 0/255])
hold on
plot(1:length(neural_id),V2(neural_id,1),'.','MarkerSize',7,...
    'MarkerFaceColor',blue,...
    'MarkerEdgeColor',blue)
ylim([-0.6,0.6])
box off
xline(0, 'k', 'LineWidth', 1);
ax = gca;        
ax.XTick = [];   
ax.XColor = 'none';
xlabel('');    
end

function plot_weight_corr(V1,V2,fit,varargin)
figure
plot(V1,V2,varargin{:})
if fit
    [beta,p_beta] = corr(V1,V2,'type','Spearman');
    hline = refline(beta,0);
    if p_beta<0.05
        hline.LineWidth = 1;
    end
    title_info = ['corr=',num2str(beta),'p=',num2str(p_beta)];
    title(title_info)
end
xlim([0,0.6])
ylim([0,0.6])
ax1.XAxisLocation='origin';
ax1.YAxisLocation='origin';
axis square
end

end

