%% plot significant time points by zn 2025/9/26
function  plotsignificant(x,y,width, height,lineColor)
    
    significant_points = x;
    y_center = y;
    if nargin < 3
        width = 5;
        height = 5; 
    end

    % 为所有点一次性计算矩形位置
    x_left = significant_points - width/2;
    x_right = significant_points + width/2;
    y_bottom = y_center - height/2;
    y_top = y_center + height/2;

    % 逐个绘制矩形
    for i = 1:length(significant_points)
        x_rect = [x_left(i), x_right(i), x_right(i), x_left(i), x_left(i)];
        y_rect = [y_bottom(i), y_bottom(i), y_top(i), y_top(i), y_bottom(i)];
        if ~isempty(x_rect)
%             plot(x_rect, y_rect, 'Color', lineColor, 'LineWidth', 2);
          fill(x_rect, y_rect, lineColor, 'FaceAlpha',  1, 'EdgeColor', lineColor, 'LineWidth', 2);
        end
    end
end