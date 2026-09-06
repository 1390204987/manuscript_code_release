% 2024/12/2 by zn find the timepoint that explained variance stop increase with
% time and return PC1 relates of data untill this time point(under mTDRanalysis)
function [seqPC1,seqtime]=getseqPC1(topEV,time,startseqtime,Xb_reshape,Bp)
[n,T] = size(Bp);
startseqindice = find(time>=startseqtime,1);
dtopEV = topEV(2:end)-topEV(1:end-1);
% find index larger than 0;
positive_indices = find(dtopEV>0);
positive_indices = positive_indices(positive_indices>startseqindice|positive_indices==startseqindice);
% find continue larger than 0 part
diff_indices = diff(positive_indices); % 计算索引差
find(diff_indices>1)
length(positive_indices)
split_points = [0 ;find(diff_indices>1);length(positive_indices)]; %分段位置
% 提取每一段连续大于 0 的索引
segments = cell(1, length(split_points) - 1);
for i = 1:length(split_points) - 1
    segments{i} = positive_indices(split_points(i)+1:split_points(i+1));
    positive_length(i,1) = length(segments{i}); 
end
which_period = find(positive_length>2,1);
if ~isempty(which_period)
    effec_period = segments{which_period};
    index_max = effec_period(end)+1;
    max_topEV = topEV(index_max);
    min_topEV = topEV(startseqindice-1);
%     if max_topEV-min_topEV<0.05% increase smaller than 5% was probably noise
%         index_max = T;
%     end
else
    index_max = T;    
end

% index_max = find(topEV == max_topEV);
if index_max~=T
    seqtime = time(index_max);
    seqdata = Bp(:,1:index_max).*Xb_reshape;
    seqdata = reshape(seqdata,n,[]);
    [seqPCs,~,~] = pca(seqdata'); %get PCs of data untill timepoint max_topEV
    seqPC1 = seqPCs(:,1);
else
    seqtime = nan;
    seqdata = Bp(:,1:index_max).*Xb_reshape;
    seqdata = reshape(seqdata,n,[]);
    [seqPCs,~,~] = pca(seqdata'); %get PCs of data untill timepoint max_topEV
    seqPC1 = seqPCs(:,1);
%     seqPC1 = nan(n,1);
end
end