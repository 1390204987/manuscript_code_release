% mTDR analysis 2024/11/3
% use mTDR for color target task analysis , rP set by youself
function mTDRLearning(datapath,savename,savepath)
%% Set Path
CrntDir = pwd;
addpath([CrntDir '/EstimatedPars'],...
    [CrntDir '/functionFiles'],...
    [CrntDir '/functionFiles/tools_kron'],...
    genpath([CrntDir '/minFunc_2012']))
if nargin==0
    datapath = 'F_STScolor_shuffleall_simult';
    savename = 'F_STScolor_shuffleall_mTDR.mat';
    savepath = './shuffle_100times/';
end
load(datapath)
P = 7;
Pb=6;
%% Parameter learning
%%
%%

%% Specify optimization parameters
% r0              = [9,10,10,12,9,4]; % F STS specify ranks.  [curr_relative_heading,curr_target_loc,curr_choice,curr_sacdir,curr_absheading]
% r0              = [9 11 11 12 10 4]; % F IPS
% r0              = [8,9,9,10,9,5]; %D STS
r0              = [7 7 7 8 8 4];% D IPS
ridgeparam      = 0;% ridge parameter for ridge regression with low-rank constraint by SVD
opts.MaxIter    = 500;
opts.Display    = 'off';
g               = 0;% ridge regression parameter for MML
rtot    = sum(r0);
%% Step 1: calculate sufficient statistics
% Start with:
% Data array Ymat (n x T x K)  where 
% n = number of neurons
% T = number of time points
% K = number of trials
% Design matrix X (K x P), where
% P = number of task variables
% Observation matrix hk (K x n) that describes which neurons were
% observed on which trials.

 %Sufficient stats for initial estimation by SVD
 % Estimation by SVD treats condition-independent component as a covariate.
 %  The consequence of this is that the design matrix has one column of all
 %  ones. i.e. X(i,:) = [var1 var2... varP 1]
[XX, XY,Yn,allstim,n,T] = MkSuffStats_BilinReg_Sims(Ymat,Xmat,hk);


% Remove condition-independent component from design matrix for MML
% estimation since we assume full rank for this component for all models
Xb = Xmat;     Xb(:,P) = [];

% Make sufficient stats for estimation by maximum marginal likelihood
ni = sum(hk,1); zi = cell(n,1); Xi = cell(n,1); Ai = zeros(Pb,Pb,n); 
Ri = zeros(Pb*T,Pb*T,n);  zzi = zeros(n,1);  Xzetai = zeros(T*Pb,n);
xbari = zeros(n,Pb); Ybar = zeros(T,n);

for ii = 1:n
    Xi{ii}          = Xb(hk(:,ii)==1,:);
    Ai(:,:,ii)      = Xi{ii}'*Xi{ii};
    zi{ii}          = squeeze(Ymat(ii,:,hk(:,ii)==1));
    zetai           = vec(zi{ii});
    Xzetai(:,ii)    = kronmult({eye(T),Xi{ii}'},zetai);
    Ri(:,:,ii)      = Xzetai(:,ii)*Xzetai(:,ii)';
    zzi(ii)         = zetai'*zetai;
    xbari(ii,:)     = mean(Xi{ii},1);
    ybar            = squeeze(sum(Ymat(ii,:,:),3))/ni(ii);
    Ybar(:,ii)      = ybar;
end

%% Step 2: initialize functions for optimization

% function for initialization by SVD
svdregress      = @(r)SVDRegress_S_Vdata(XX,XY,Yn,allstim,T,n,r,ridgeparam,opts);

% function for direct optimization of marginal likelihood
MMLEParfun      = @(lamb0,s0,bhat0,r)Estpars_CoordAscent_lambi_S_b(lamb0,s0,bhat0,r,Ai,Xi,zi,ni,xbari,Ybar,Xzetai);

% function for optimization by EM
EMregressfun    = @(r,pars0)ECMEtdr('converge',1e-0,pars0,Ai,Xi,r,ni,zi,xbari,Ybar,Xzetai);
EMparsfun       = @(r)ECMEregress_wrapper(r,svdregress,EMregressfun,Ybar);

%% Step 3: run optimization and parse parameter vector

% function implementing end-to-end optimization with SVD initialization
MMLE_EMregressfun   = @(r)MMLE_CoordAscentWrapper(EMparsfun,MMLEParfun,r,n,T);
parhist            = MMLE_EMregressfun(r0);
histfileMMLE = 'EstimatedPars/LearnDemoMMLE';
save(histfileMMLE,'parhist')

lambhat = parhist(1:n);% estimated noise precisions
Bhat0   = reshape(parhist(n+rtot*T+1:end),T,n);% estimated condition-independent components
[~,~,Xzetai_opt] = ECMEsuffstat(zi,Xi,Bhat0);% recalculated a sufficient stat
[Bhat, r,What,Shat,lambhat] = MakeBhat_data(histfileMMLE,Ai,Xzetai_opt,r0,[]);% estimated low-rank regression parameters

save([savepath savename],'Bhat', 'r','What','Shat','lambhat','Bhat0')


%%  Seq PCA
%{
close all
time = t_centers{:};
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
     color_{6,1} = [0,0,1;0,0,0.9;0,0,0.8;0,0,0.7;0,0,0.6]; %coherence
     color_{7,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
     linesize = 2;
     color = color_{p,1};
     variable_name={'heading';'color';'choice';'saccade';'last_sac';'coherence'};
    topEV{1,1} = gettopEV(Xb_reshape,plotB,startlen,startpoint);   
    stupid_time = find(time>-1000&time<800);
    startseqtime = time(find(topEV{1,1}==min(topEV{1,1}(stupid_time))));
    [seqPC(:,1),seqtime1]=getseqPC1(topEV{1,1},time,startseqtime,Xb_reshape,plotB);

    subperiod1_B = plotB-seqPC(:,1)*seqPC(:,1)'*plotB;    
    topEV{2,1} = gettopEV(Xb_reshape,subperiod1_B,startlen,startpoint);
    startseqtime = seqtime1;
    [seqPC(:,2),seqtime2]=getseqPC1(topEV{2,1},time,startseqtime,Xb_reshape,subperiod1_B);
    
    subperiod2_B = subperiod1_B-seqPC(:,2)*seqPC(:,2)'*subperiod1_B;
    topEV{3,1} = gettopEV(Xb_reshape,subperiod2_B,startlen,startpoint);
    startseqtime = seqtime2;
    [seqPC(:,3),seqtime3]=getseqPC1(topEV{3,1},time,startseqtime,Xb_reshape,subperiod2_B);
    
    %plot psth projection on seqPC(:,1)
%     [psth_heading,unique_cons,unique_cohe] = get_PSTH(Xmat,Ymat,1,hk);
    set(figure,'position',[60 91 900 300]);
    plotpsth1D(time,psth,unique_cons,unique_cohe,seqPC(:,1),color,linesize)   
    ylims = ylim;
    plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
    plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
    ax = gca;
%     ax.YTick = [];  
%     yticklabels([]);
    ylabel([variable_name{p},'axis1'])
    SetFigure
   %plot divergence projection on seqPC(:,1)
    set(figure,'position',[60 91 900 300]);
    seqPC1_diverge(:,:,p) = plot_divergence(time,psth,unique_cons,unique_cohe,seqPC(:,1),color,linesize);  
    ylims = ylim;
    plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
    plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
    ax = gca;
    ylabel([variable_name{p},'axis1'])
    SetFigure
    
    %plot psth projection on seqPC(:,2)
    set(figure,'position',[60 91 900 300]);
    plotpsth1D(time,psth,unique_cons,unique_cohe,seqPC(:,2),color,linesize)
    ylims = ylim;
    plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
    plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
    ax = gca;
%     ax.YTick = [];  
%     yticklabels([]);
    ylabel([variable_name{p},'axis2'])
    SetFigure
    %plot divergence projection on seqPC(:,2)
    set(figure,'position',[60 91 900 300]);
    seqPC2_diverge(:,:,p)=plot_divergence(time,psth,unique_cons,unique_cohe,seqPC(:,2),color,linesize);
    ylims = ylim;
    plot([seqtime1,seqtime1],ylims,'k','LineWidth',2)
    plot([seqtime2,seqtime2],ylims,'k','LineWidth',2)
    ax = gca;
    ylabel([variable_name{p},'axis2'])
    SetFigure
    
    %plot psth projection on seqPC(:,3)
    figure
    plotpsth1D(time,psth,unique_cons,unique_cohe,seqPC(:,3),color,linesize)
    ylabel([variable_name{p},'axis3'])
    SetFigure
    %plot divergence projection on seqPC(:,3)
    figure
    seqPC3_diverge(:,:,p)=plot_divergence(time,psth,unique_cons,unique_cohe,seqPC(:,3),color,linesize);
    ylabel([variable_name{p},'axis3'])
    SetFigure
    
    
    variablePC{p} = [seqPC(:,1),seqPC(:,2),seqPC(:,3)];
    
    % plot top EV 
    set(figure,'position',[60 91 900 300]);
    axis off
    title_infor = ['variable P  = ',num2str(p)];
    title(title_infor)
    h_subplot = tight_subplot(1,3,[0.01 0.02],[0.1 0.1],[0.08 0.03],[],[]);    
    for i = 1:3
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
    SetFigure
    % plot B projection

    set(figure,'position',[60 91 900 300]);
    axis off
    title_infor = ['variable P  = ',num2str(p)];
    title(title_infor)
    h_subplot = tight_subplot(1,3,[0.01 0.02],[0.1 0.1],[0.08 0.03],[],[]);    
    for i = 1:3
        set(gcf,'CurrentAxes',h_subplot(i))
        hold on
        B_project = seqPC(:,i)'*plotB;
        plot(time,B_project)
        ylims = ylim;
        plot([seqtime1,seqtime1],ylims)
        plot([seqtime2,seqtime2],ylims)
        plot([seqtime3,seqtime3],ylims)
    end
    SetFigure
    figure
    plot2D(time,psth,unique_cons,unique_cohe,seqPC(:,1),seqPC(:,2),color,linesize)    
    ax = gca;
    ax.YTick = [];
    yticklabels([]);
    ax.XTick = [];
    xticklabels([]);
    xlabel('choice axis 1')
    ylabel('choice axis 2')
    SetFigure
    
    figure
    plot3D(time,psth,unique_cons,unique_cohe,seqPC(:,1),seqPC(:,2),seqPC(:,3),color,linesize)    
    SetFigure
 end
%plot each variable divergence together
 figure
 plotdivergeall(time,seqPC1_diverge) 
 SetFigure
%  save(savename,'variablePC','cellID_pick')
a=1
%%
% load('D_variablePC_ips.mat')
% close all
index = randperm(length(variablePC{1}(:,1))); 
V1 = variablePC{1}(:,1);
V2 = variablePC{3}(:,1);
[~,neural_id] = sort(V1);
set(figure,'position',[60 91 600 100]);
SetFigure
plot_weight(V1,V2,neural_id);

plot_weight_corr(V1,V2,1,'o','MarkerEdgeColor','k','MarkerFaceColor','k')

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
%}
end
