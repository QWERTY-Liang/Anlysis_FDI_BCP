function EMG_processed = TL6_preprocessEMG(EMGdata, fs)
% Author: Liang Tong
% Date: 10/9/2024     
% 输入：
    % EMGdata: 原始EMG数据
    % Fs: 采样频率
    %
    % 输出：
    % EMG_processed: 经过处理的EMG数据

      % 输入：
    % EMGdata: 1×3 cell array， 每个cell包含两列数据（左、右）
    % Fs: 采样频率
    %
    % 输出：
    % EMG_processed: 1×3 cell array，处理后的数据

    EMG_processed = cell(1, 3);  % 初始化输出 cell array
    
%% 0. 将第三个 cell 中的 NaN 值替换为 0
    EMGdata{3}(isnan(EMGdata{3})) = 0;

    %% 1. 50Hz 陷波滤波器（Notch Filter）和 50-150Hz 带通滤波器（Bandpass Filter）
    % wo = 50 / (Fs/2);  % 50 Hz 标准化频率
    % bw = wo / 35;      % 带宽，较小的值使得滤波器更窄
    % [b_notch, a_notch] = iirnotch(wo, bw);
    % 
    % % 定义带通滤波器参数
    % low_cutoff = 50;  % 带通滤波器低截止频率
    % high_cutoff = 150; % 带通滤波器高截止频率
    % [b_bp, a_bp] = butter(4, [low_cutoff, high_cutoff] / (Fs/2), 'bandpass');
     % 设计FIR带通滤波器
         low_cutoff = 50; % 低截止频率
    high_cutoff = 150; % 高截止频率
    filter_order = 4; % 滤波器阶数
    bp_filter = designfilt('bandpassfir', 'FilterOrder', filter_order, ...
                           'CutoffFrequency1', low_cutoff, 'CutoffFrequency2', high_cutoff, ...
                           'SampleRate', fs);



    f0_50Hz = 50; % Frequency to notch (in Hz)
Q = 35; % Quality factor (higher value means a narrower notch)

% Design the 50 Hz notch filter
notchFilter_50Hz = designfilt('bandstopiir', ...
                              'FilterOrder', 2, ...
                              'HalfPowerFrequency1', f0_50Hz*(sqrt(2)/(2*Q)), ...
                              'HalfPowerFrequency2', f0_50Hz*(sqrt(2)/(2*Q)) * (Q+1)/Q, ...
                              'DesignMethod', 'butter', ...
                              'SampleRate', fs);

    % 对所有 cell 中的数据进行滤波
    for i = 1:3
        for j = 1:2  % 处理左、右两列
            data = EMGdata{i}(:, j);
            % 50Hz 陷波滤波
            data_notched = filtfilt(notchFilter_50Hz, data);
            % 50-150Hz 带通滤波
            data_bandpassed = filtfilt(notchFilter_50Hz, data_notched);
            EMG_processed{i}(:, j) = data_bandpassed;  % 储存滤波后的数据
        end
    end

    %% 2. 均值中心化 (Mean Centering)
    % 使用第一个cell的均值进行中心化
    baseline_mean_left = mean(EMG_processed{1}(:, 1));  % 第一个cell左列均值
    baseline_mean_right = mean(EMG_processed{1}(:, 2)); % 第一个cell右列均值

    for i = 1:3
        EMG_processed{i}(:, 1) = EMG_processed{i}(:, 1) - baseline_mean_left;
        EMG_processed{i}(:, 2) = EMG_processed{i}(:, 2) - baseline_mean_right;
    end

    %% 3. 整流 (Rectification)
    for i = 1:3
        EMG_processed{i}(:, 1) = abs(EMG_processed{i}(:, 1));  % 左列整流
        EMG_processed{i}(:, 2) = abs(EMG_processed{i}(:, 2));  % 右列整流
    end

    %% 4. Z-Score 标准化 (Z-Score Normalization)
    % % 使用第一个cell的数据进行标准化
    baseline_std_left = std(EMG_processed{1}(:, 1));  % 第一个cell左列标准差
    baseline_std_right = std(EMG_processed{1}(:, 2)); % 第一个cell右列标准差

    for i = 1:3
        EMG_processed{i}(:, 1) = EMG_processed{i}(:, 1) / baseline_std_left;  % 左列标准化
        EMG_processed{i}(:, 2) = EMG_processed{i}(:, 2) / baseline_std_right; % 右列标准化
    end

   %% 5.1  Low pass filter （smoothoing method 1） 二选一
    % 使用 designfilt 设计 15 Hz 低通滤波器 （越小delay越高，越小形状越清晰）
    % cutoff=15;
    % low_pass_filter = designfilt('lowpassiir', 'FilterOrder', 4, ...
    %                              'HalfPowerFrequency', cutoff/(fs/2), ...
    %                              'DesignMethod', 'butter');
    % 
    % % 对所有 cell 应用低通滤波
    % for i = 1:3
    %     for j = 1:2  % 处理左、右两列
    %         EMG_processed{i}(:, j) = filtfilt(low_pass_filter, EMG_processed{i}(:, j));
    %     end
    % end

    %% 5.2 RMS smoothing (smoothing method 2) 二选一
    % window_size = round(0.1 * fs);  % 100 ms 窗口
    % for i = 1:3
    %     for j = 1:2  % 处理左、右两列
    %         EMG_processed{i}(:, j) = rms_smooth(EMG_processed{i}(:, j), window_size);
    %     end
    % end
end

%% RMS 平滑函数
function data_rms = rms_smooth(data, window_size)
    % 定义滑动窗口的 RMS 计算
    data_rms = sqrt(movmean(data.^2, window_size));
end
