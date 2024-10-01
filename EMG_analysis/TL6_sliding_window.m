function results = TL6_sliding_window(datafull, window_size, step_size)
m=2;r=0.4;   
% 获取datafull的长度
    data_length = length(datafull);
    
    % 计算滑动窗口的数量
    num_windows = floor((data_length - window_size) / step_size) + 1;
    
    % 初始化结果数组
   results = zeros(num_windows, 1);
     % 生成汉宁窗
    taper_window = hann(window_size);

    % 滑动窗口处理
    for i = 1:num_windows
        % 窗口的起始和结束位置
        start_idx = (i - 1) * step_size + 1;
        end_idx = start_idx + window_size - 1;
        
        % 提取窗口数据
        window_data = datafull(start_idx:end_idx);
         % 应用汉宁窗
        tapered_window_data = window_data .* taper_window;
        
        % 调用自己的程序处理窗口数据并存储结果
        results(i) = TL6_FuzzyEn( tapered_window_data,m,r);
    end
end

