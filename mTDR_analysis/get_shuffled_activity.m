% get shuffleed data projetioned activity
function [shu_activity] = get_shuffled_activity(monkey,brain_area,p)
    if nargin==0
        monkey = ['D'];
        brain_area = ['STS'];
        p=3;
    end
        vec_path = ['./shuffle_20times/',monkey,'_',brain_area,'/'];
        data_path = ['./shuffle_20times/',monkey,'_',brain_area,'/'];
        color_{1,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
        color_{2,1} = [0,191,255;155,48,255]/255; %color
        color_{3,1} = [139,10,80;72,61,139]/255; %choice
        color_{4,1} = [0,0,128;205,133,0]/255; %saccade
        color_{5,1} = [0,0,139;205,133,0]/255; %last saccade
        color_{6,1} = [0,0,1;0,0,0.9;0,0,0.8;0,0,0.7;0,0,0.6]; %coherence
        color_{7,1} = [205,0,0;0,0,0;0,100,0]/255; %heading
        color = color_{p,1};
        linesize = 2;

        for ishuffle = 1:10
            vec_name = [vec_path,monkey,'_colorvariable_shuffleall_all_',num2str(ishuffle),'_',brain_area];
            data_name = [data_path,monkey,'_',brain_area,'color_shuffleall_all_20ms_',num2str(ishuffle),'.mat'];
            shuffle_vec = load(vec_name); 
            shuffle_data = load(data_name);
            
            shu_time = shuffle_data.t_centers{:}; 
            if p == 1
                [shu_psth,heading_cons,choice_cons,unique_cohe] = get_2v_PSTH(shuffle_data.Xmat,shuffle_data.Ymat,shuffle_data.hk);
                unique_cons = heading_cons;
                control_cons = choice_cons;
                shu_psth = permute(shu_psth,[1,2,4,3]);
            elseif p == 3
                [shu_psth,heading_cons,choice_cons,unique_cohe] = get_2v_PSTH(shuffle_data.Xmat,shuffle_data.Ymat,shuffle_data.hk);
                unique_cons = choice_cons;
                control_cons = heading_cons;
            else
                [shu_psth,~,~,unique_cons,unique_cohe] = get_PSTH(shuffle_data.Xmat,shuffle_data.Ymat,p,shuffle_data.hk);
                control_cons = [];
            end
            icon1 = find(unique_cons==-1);
            icon2 = find(unique_cons==1);
            coherence_0=0; 
            
            [mean_shu_divergepsth_seqPC1]=getdivergence_shuffle(shu_time,shu_psth,...
            unique_cons,control_cons,unique_cohe,shuffle_vec.variablePC{1,p}(:,1),coherence_0,color,linesize);
            temp = mean_shu_divergepsth_seqPC1;
            abs_temp = abs(temp);
            maxindex = find(abs_temp == max(abs_temp));
            if temp(maxindex)>0
                temp = temp;
            else
                temp = -temp;
            end
%             shu_activity(ishuffle,:) = mean_shu_divergepsth_seqPC1;       
            shu_activity(ishuffle,:) = temp;       
        end
        mean_shu_ac = mean(shu_activity);
        sd_shu_ac = std(shu_activity);
        figure()
        shadedErrorBar(shu_time,mean_shu_ac,sd_shu_ac,...
        'lineprops',{'Color',color(1,:),'LineStyle','--'},'transparent',1);
        ylim([-0.1,0.1])
end