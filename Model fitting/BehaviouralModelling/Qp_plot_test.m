% 生成模拟反应时间数据 (单位：秒)
RT_correct = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0]; % 正确反应
RT_error   = [0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0];  % 错误反应

% 生成相应的响应概率
p_correct = linspace(0, 1, length(RT_correct)); % 正确反应概率
p_error   = linspace(0, 1, length(RT_error));   % 错误反应概率

% 计算 RT 分位数（例如 10%、30%、50%、70%、90%）
quantiles = [0.1, 0.3, 0.5, 0.7, 0.9];

RT_correct_quantiles = quantile(RT_correct, quantiles);
RT_error_quantiles   = quantile(RT_error, quantiles);

figure; hold on;

% 绘制正确反应的 RT 分位数
plot(p_correct, RT_correct_quantiles, '-o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'DisplayName', 'Correct');

% 绘制错误反应的 RT 分位数
plot(1 - p_error, RT_error_quantiles, '-s', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'DisplayName', 'Error');

% 设置坐标轴和图例
xlabel('Response Probability');
ylabel('RT Quantile (s)');
title('Quantile-Probability Plot');
legend('show');
grid on;