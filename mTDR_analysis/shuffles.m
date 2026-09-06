% generate 100times shuffle 
shuffle_number = 10;
for ishuffle = 1: shuffle_number
    FolderPath_colorHD =   {'Z:\Data\Tempo\Batch\fingerpostDrogon_LA\DrogoncolorcoarseSTS_LA\';};
% FolderPath_colorHD = {'Z:\Data\Tempo\Batch\fingerpostFriday_LA\FridaycolorIPS_LA\';'Z:\Data\Tempo\Batch\fingerpostFriday_NP\FridaycolorIPS_NP\'};
    
    savepath = './shuffle_20times/D_STS/';
    time_bin = 100; %ms
    savename = ['D_STScolor_shuffleall_all_',num2str(time_bin),'ms_',num2str(ishuffle)];
    shuffle_trial = 1;
    randseed = ishuffle;
    data_rerange(FolderPath_colorHD, time_bin, savepath,savename, shuffle_trial,randseed)
    
    time_bin = 20; %ms
    savename = ['D_STScolor_shuffleall_all_',num2str(time_bin),'ms_',num2str(ishuffle)];
    data_rerange(FolderPath_colorHD, time_bin, savepath,savename, shuffle_trial,randseed)
end