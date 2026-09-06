function plot_roc_cor(var1,var2,fit,text_y,varargin)
% varargin 1 plot point parameter
% varargin 2 plot line parameter

% axis equal
xlim([0,1]);
ylim([0,1]);

if ~isempty(varargin{2})
    plot(var1,var2,varargin{1}{:},'MarkerEdgeColor',varargin{2}{:})
else
    plot(var1,var2,varargin{1}{:})
end
if fit
    valid_indices = ~isnan(var1) & ~isnan(var2); 
    [beta ,p_beta] = corr(var1(valid_indices),var2(valid_indices),'type','Pearson');
    intercept = 0.5-0.5*beta;
%     hline = refline(beta,intercept);

%     x_data = var1(valid_indices);
%     y_data = var2(valid_indices);
%     X = [ones(length(x_data), 1), x_data];
%     % 线性回归
%     [b, bint, r, rint, stats] = regress(y_data, X);
%     beta = b(2);        % 斜率
%     intercept = b(1);    % 截距
%     p_beta = stats(3);  % p 值
        
%     hline = refline(beta, intercept);
    if p_beta>0.05|isnan(p_beta)
        hline.LineWidth = 0.5;
        hline.LineStyle = ':';
    else
        hline.LineWidth = 2;
        hline.LineStyle = '-';
        if ~isempty(varargin{2})
            hline.Color = varargin{2}{:};
        end
    end    
    text(0.6,text_y+0.1,['corr=', sprintf('%.3g', beta)],'Color',varargin{2}{:})
    text(0.6,text_y,['p=',  sprintf('%.3g', p_beta)],'Color',varargin{2}{:})
end
plot([0.5 0.5],[0,1],'--k')
plot([0.0 1],[0.5,0.5],'--k')
plot([0.1 0.9],[0.9,0.1],'-k')
plot([0.1 0.9],[0.1,0.9],'-k')
end