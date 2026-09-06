%% 2025/2/8 get each psth of color target HD task(under single trial) for unique conditon and unique coherence
function [psth,psth_correct,psth_error,unique_con,unique_coherence] = get_PSTH(Xmat,Ymat,variable_id,hk)
var_info = Xmat(:,variable_id);

if size(Xmat,2)==7 % color HD task
    if variable_id ==1 % heading
        unique_con = [-1,0,1];
    elseif variable_id == 6 % coherence
        unique_con = [1];
    else
        unique_con = unique(var_info);
    end
    coherence_info = Xmat(:,6);
    unique_coherence =  unique(coherence_info);
    outcome_info = sign(Xmat(:,1)).*sign(Xmat(:,3));
elseif size(Xmat,2)==5 % HD task
    if variable_id ==1 % heading
        unique_con = [-1,0,1];
    elseif variable_id == 4 % coherence
        unique_con = [1];
    else
        unique_con = unique(var_info);
    end
    coherence_info = Xmat(:,4);
    unique_coherence =  unique(coherence_info);
    outcome_info = sign(Xmat(:,1)).*sign(Xmat(:,2));
end
unique_outcome = [-1,1];
% normallize psth
% Z_Ymat= (Ymat-mean(Ymat,3))./std(Ymat,0,3); can't use this as some neuron
% doesn't have activity under some trial.
[~,n] = size(hk);
[~,T,~] = size(Ymat);
Z_Ymat = zeros(size(Ymat)); % ³õÊ¼»¯Îª 0
for ineuro = 1:n
    std_val = nanstd(Ymat(ineuro,:,hk(:,ineuro)), 0, 3);
    non_zero_std = std_val ~= 0; % for some timepints firing_rate==0
    y_neuron = Ymat(ineuro,:,hk(:,ineuro));
        Z_Ymat(ineuro, non_zero_std,:) = (Ymat(ineuro, non_zero_std,:) - nanmean(Ymat(ineuro, non_zero_std, hk(:,ineuro)), 3))./ nanstd(y_neuron,0,3);
%     Z_Ymat(ineuro, non_zero_std,:) = (Ymat(ineuro, non_zero_std,:) - nanmean(Ymat(ineuro, non_zero_std, hk(:,ineuro)), 3))./ nanstd(y_neuron(:));
%             Z_Ymat(ineuro, non_zero_std,:) = (Ymat(ineuro, non_zero_std,:) - nanmean(Ymat(ineuro, non_zero_std, hk(:,ineuro)), 3))./max(max(Ymat(ineuro,:, hk(:,ineuro))));
%         end
     for icohes = 1:length(unique_coherence)
        select_coh = coherence_info == unique_coherence(icohes);
        for icon = 1:length(unique_con)
            
            if variable_id==1
                select = sign(var_info)==unique_con(icon);
                if unique_coherence(icohes)~=0
                    select_correct = select&outcome_info==1;
                    select_error = select&outcome_info==-1;
                else
                    select_correct = select;
                    select_error = select;
                end
                %                 select = outcome_info==unique_outcome(icon)&sign(var_info)==1;% conditioned according to outcome
            elseif variable_id==3
                %                 select = outcome_info==unique_outcome(icon);% conditioned according to outcome
                select = var_info==unique_con(icon);% conditioned according to choice
                if unique_coherence(icohes)~=0
                    select_correct = select&outcome_info==1;
                    select_error = select&outcome_info==-1;
                else
                    select_correct = select;
                    select_error = select;
                end       
            elseif variable_id==6
                select = ones(length(var_info),1);
                if unique_coherence(icohes)~=0
                    select_correct = select&outcome_info==1;
                    select_error = select&outcome_info==-1;
                else
                    select_correct = select;
                    select_error = select;
                end
            else
                select = var_info==unique_con(icon);
                if unique_coherence(icohes)~=0
                    select_correct = select&outcome_info==1;
                    select_error = select&outcome_info==-1;
                else
                    select_correct = select;
                    select_error = select;
                end
            end
            %             select_acti = squeeze(nanmean(Z_Ymat(ineuro,:,select&hk(:,ineuro)&select_coh),3));

            select_acti = reshape(Z_Ymat(ineuro,:,select&hk(:,ineuro)&select_coh),T,[]);
%             select_acti = reshape(Ymat(ineuro,:,select&hk(:,ineuro)&select_coh),39,[]);
            correct_select_acti = reshape(Z_Ymat(ineuro,:,select_correct&hk(:,ineuro)&select_coh),T,[]);
%             correct_select_acti = reshape(Ymat(ineuro,:,select_correct&hk(:,ineuro)&select_coh),39,[]);
            if sum(select_error&hk(:,ineuro)&select_coh)~=0
                error_select_acti = reshape(Z_Ymat(ineuro,:,select_error&hk(:,ineuro)&select_coh),T,[]);
%                 error_select_acti = reshape(Ymat(ineuro,:,select_error&hk(:,ineuro)&select_coh),39,[]);
            else
                error_select_acti = nan(T,1);
            end
            psth{ineuro,icohes,icon} = select_acti;
            psth_correct{ineuro,icohes,icon} = correct_select_acti;
            psth_error{ineuro,icohes,icon} = error_select_acti;
        end
    end
end
end