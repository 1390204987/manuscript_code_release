% by zn 2025/4/5 analysis axis correlation get by mTDR and seqPCA
clear; close all;
%%
Vec_mTDR = load('D_IPScolor_all_mTDR.mat');
Vec_seqPCA = load('D_colorvariable_all_IPS_redit.mat');
% Vec_seqPCA = load('D_colorvariable_all_IPS.mat');
heading_vec = Vec_seqPCA.variablePC{1};
% heading_dim =  ~any(isnan(heading_vec), 1); % 找出不含 NaN 的列（逻辑索引）
% heading_dim = Vec_seqPCA.real_dim{1}(:);
heading_dim = 1;
heading_vec_name = {'heading_vec1','heading_vec2','heading_vec3'};
color_vec = Vec_seqPCA.variablePC{2};
% color_dim =  ~any(isnan(color_vec), 1); % 找出不含 NaN 的列（逻辑索引）
% color_dim = Vec_seqPCA.real_dim{2}(:);
color_dim = 1;
color_vec_name = {'color_vec1','color_vec2','color_vec3'};
choice_vec = Vec_seqPCA.variablePC{3};
% choice_dim =  ~any(isnan(choice_vec), 1); % 找出不含 NaN 的列（逻辑索引）
choice_dim = Vec_seqPCA.real_dim{3}(:);
choice_vec_name = {'choice_vec1','choice_vec2','choice_vec3'};
sac_vec = Vec_seqPCA.variablePC{4};
% sac_dim =  ~any(isnan(sac_vec), 1); % 找出不含 NaN 的列（逻辑索引）
sac_dim =  Vec_seqPCA.real_dim{4}(:);

sac_vec_name = {'sac_vec1','sac_vec2','sac_vec3'};

allvec = [heading_vec(:,heading_dim),color_vec(:,color_dim),choice_vec(:,choice_dim),sac_vec(:,sac_dim)];
allvec_name = {heading_vec_name{heading_dim},color_vec_name{color_dim},choice_vec_name{choice_dim},sac_vec_name{sac_dim}};
[corr_mat,p_mat] = corr(allvec,'Type','pearson');
vecdot_mat = allvec'*allvec;
norm_vec = sqrt(sum(allvec.^2));
norm_mat = norm_vec'*norm_vec;
cos_angle = vecdot_mat./norm_mat;
angle_rad = acos(cos_angle);
 angle_deg = real(rad2deg(angle_rad));
% angle_deg = min(real(rad2deg(angle_rad)),180-real(rad2deg(angle_rad)));


% Shuffle test for pairwise axis angles, 0-180 deg
nShuffle = 1000;

[p_angle, shuffle_angle, shuffle_angle_mean, shuffle_angle_CI, real_angle_shuffle] = ...
    angle_shuffle_test(allvec, nShuffle);

% Heatmap of shuffle-based p values
figure()

p_plot = p_angle;
p_plot(p_plot > 0.05) = NaN;

h3 = heatmap(allvec_name, allvec_name, p_plot);
h3.ColorLimits = [0, 0.05];
h3.Colormap = parula;
h3.MissingDataColor = 'white';

title('Shuffle-based p value of axis alignment');

SetFigure



% define colormap
% 从天蓝到浅蓝再到白色
key_colors = [0.53 0.81 1;   % 天蓝色
              0.75 0.90 1;   % 浅蓝色
              1    1    1];  % 白色  
n = 256;
custom_cmap = interp1(linspace(0, 1, size(key_colors, 1)), key_colors, linspace(0, 1, n));

% heat map of corr_mat
figure()
p_mat(p_mat > 0.05) = NaN; % 将 >0.05 的值设为 NaN
h1 = heatmap(allvec_name, allvec_name, p_mat);
h1.ColorLimits = [0, 0.05]; % 强制颜色映射到 0-0.05
h1.Colormap = parula; % 使用 parula 颜色方案（或其他如 jet, hot 等）
h1.MissingDataColor = 'white'; % NaN 显示为白色
SetFigure
% heat map of axis angle
figure()
h2 = heatmap(allvec_name, allvec_name, angle_deg);
h2.ColorLimits = [0, 90]; % 强制颜色映射到 0-0.05
h2.Colormap = custom_cmap; % 使用 parula 颜色方案（或其他如 jet, hot 等）
h2.MissingDataColor = 'white'; % NaN 显示为白色
SetFigure

% heading axis compare
% figure()
% plot_pearson_corr(heading_vec(:,1),choice_vec(:,1),1,'ob')
% SetFigure
% figure()
% plot_pearson_corr(abs(heading_vec(:,1)),abs(sac_vec(:,1)),1,'ok')
% SetFigure
figure()
% plot_pearson_corr(choice_vec(:,1),sac_vec(:,1),1,'^k')
plot_pearson_corr(sac_vec(:,1),choice_vec(:,1),1,'^k')
SetFigure
figure()
plot_pearson_corr(heading_vec(:,1),choice_vec(:,1),1,'sk')
SetFigure
% 
% figure()
% plot_pearson_corr(abs(color_vec(:,1)),abs(choice_vec(:,1)),1,'sk')
% SetFigure
