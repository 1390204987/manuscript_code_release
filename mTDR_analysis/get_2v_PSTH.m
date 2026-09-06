%% 2025/4/12 get each psth of color target HD task(under single trial) for 2 conditon and unique coherence
% this is specific for heading and choice divergence calculation
function [psth,unique_heading_con,unique_choice_con,unique_coherence] = get_2v_PSTH(Xmat,Ymat,hk)
if size(Xmat,2)==7 % color HD task
    heading_info = Xmat(:,1);
    choice_info = Xmat(:,3);
    unique_heading_con = [-1,0,1];
    unique_choice_con = [-1,1];
    coherence_info = Xmat(:,6);
elseif size(Xmat,2)==5 %HD task
    heading_info = Xmat(:,1);
    choice_info = Xmat(:,2);
    unique_heading_con = [-1,0,1];
    unique_choice_con = [-1,1];
    coherence_info = Xmat(:,4);
end
    
    unique_coherence =  unique(coherence_info);
    [~,n] = size(hk);
    [~,T,~] = size(Ymat);
    
    for ineuro = 1:n
        std_val = nanstd(Ymat(ineuro,:,hk(:,ineuro)), 0, 3);
        non_zero_std = std_val ~= 0; % for some timepints firing_rate==0
        Z_Ymat = zeros(size(Ymat(ineuro,:,:))); % ³õÊ¼»¯Îª 0
        y_neuron = Ymat(ineuro,:,hk(:,ineuro));
        Z_Ymat(ineuro, non_zero_std,:) = (Ymat(ineuro, non_zero_std,:) - nanmean(Ymat(ineuro, non_zero_std, hk(:,ineuro)), 3))./ nanstd(y_neuron(:));
        for icohes = 1:length(unique_coherence)
            select_coh = coherence_info == unique_coherence(icohes);
            for icon_heading = 1:length(unique_heading_con)
                for icon_choice = 1:length(unique_choice_con)
                    select_heading = sign(heading_info)==unique_heading_con(icon_heading);
                    select_choice = sign(choice_info)==unique_choice_con(icon_choice);
                    select_acti = reshape(Z_Ymat(ineuro,:,select_heading&select_choice&hk(:,ineuro)&select_coh),T,[]);
                    psth{ineuro,icohes,icon_heading,icon_choice} = select_acti;
                end
            end
        end
   end