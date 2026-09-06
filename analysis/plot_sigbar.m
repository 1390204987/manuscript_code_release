%% plot significant bar
function plot_sigbar(x,y,step_size,color)
    % 找出连续的分组
    diff_x = diff(x);  % 计算相邻x的差值
    group_breaks = find(diff_x > step_size);  % 找出不连续的位置（差值>1）
    group_starts = [1, group_breaks + 1];  % 每组的起始索引
    group_ends = [group_breaks, length(x)];  % 每组的结束索引
    
        % 对每一组连续的点绘制线段
    for i = 1:length(group_starts)
        idx = group_starts(i):group_ends(i);  % 当前组的索引
        plot(x(idx), y(idx), ...
            '-', 'Color', color, 'LineWidth', 4);
    end
end