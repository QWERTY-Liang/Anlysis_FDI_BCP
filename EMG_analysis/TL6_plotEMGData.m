function TL6_plotEMGData(EMGdata, Fs)
    % 输入：
    % EMGdata: 1×3 cell array， 每个cell包含两列数据（左、右）
    % Fs: 采样频率
    
    % 获取三个cell的数据长度
    len1 = size(EMGdata{1}, 1);
    len2 = size(EMGdata{2}, 1);
    len3 = size(EMGdata{3}, 1);

    % 构造时间轴
    time1 = (0:len1-1) / Fs;
    time2 = (len1:len1+len2-1) / Fs;
    time3 = (len1+len2:len1+len2+len3-1) / Fs;

    % 将时间和数据拼接在一起
    time_total = [time1, time2, time3];
    left_data = [EMGdata{1}(:, 1); EMGdata{2}(:, 1); EMGdata{3}(:, 1)];  % 左手
    right_data = [EMGdata{1}(:, 2); EMGdata{2}(:, 2); EMGdata{3}(:, 2)]; % 右手

    % 创建图形窗口
    figure;

    % 绘制左手和右手数据在同一图中
    plot(time_total, left_data, 'r', 'DisplayName', 'Left Hand');
    hold on;
    plot(time_total, right_data, 'b', 'DisplayName', 'Right Hand');
    hold on;

    % 设置图表标题和标签
    title('Left and Right Hand EMG Data');
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % 标记不同阶段的分界线
    xline(time1(end), '--k', 'End of Baseline');
    xline(time2(end), '--k', 'End of Trial');

    % 确保左右手数据具有相同的 y 轴限制
    ylim([min([left_data; right_data]), max([left_data; right_data])]);

    % 添加图例
    legend;
    
    hold off;
end
