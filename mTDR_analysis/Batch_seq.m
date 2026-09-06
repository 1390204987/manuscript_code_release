% batch to calculate seq PCA
shuffle_num = 10;
for i = 1:shuffle_num
    data_path = ['./shuffle_20times/D_STS/'];
    TDR_name = [data_path,'D_STScolor_shuffleall_all_',num2str(i),'_mTDR.mat'];
    data_name = [data_path,'D_STScolor_shuffleall_all_100ms_',num2str(i),'.mat'];
    plot_fig = 0;
    save_path = ['./shuffle_20times/D_STS/'];
    save_name = ['D_colorvariable_shuffleall_all_',num2str(i),'_STS.mat';];
    Seq_PCA(TDR_name, data_name,plot_fig,save_path,save_name)
end