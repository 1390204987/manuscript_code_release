% get axis for each shuffled data
shuffle_num = 10;
for i = 1:shuffle_num
    data_name = ['D_IPScolor_shuffleall_all_100ms_',num2str(i),'.mat'];
    data_path = ['./shuffle_20times/D_IPS/',data_name];
    save_name = ['D_IPScolor_shuffleall_all_',num2str(i),'_mTDR.mat'];
    save_path = ['./shuffle_20times/D_IPS/'];
    mTDRLearning(data_path,save_name,save_path)
end