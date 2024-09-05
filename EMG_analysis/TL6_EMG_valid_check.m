%%  6.1 EMG analysis
% Author: Liang Tong
% Date: 5/9/2024 

%最后确认更新有效性，手动检查步骤
%% Toolbox requirements: 
clc
clear all
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2'));% current folder

%% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);


%% 1. Load rl_BB pre-cue-188
% col 1: subject num
% col 2: block num
% col 3: contrast (high or low contrast)
% col 4: muscle (1=FDI, 2=BCP)
% col 5: trial outcome 345 is useless(1=correct, 2=error, 4=no response, 3 too early, 6 slow, 5 wrong muscle)
% col 6: RT in sec
% col 7: evshowtime =[]   % 0.2sec additional time after response
% col 8: participant response, 1 = left, 2 = right 0= on response
% col 9: correct response, 1 = left, 2 = right
% col 10: first tilt used for SSVEP

%new EMG onsite data created
%col 11 是更新后的 EMG onsite time
%col 12 是与原来时间的提前量
%col 13 是左手brust计数
%col 14 是右手brust计数
%col 15 是是否更新 1为更新，0，-1为无更新或太长
addpath(genpath(exp.finalpath));

% load TL_ALL_include_EMG_updated.mat% EMG missing BCP totEMG
% 
% load TL_ALL_include_EMG% in EMGfolder

load([exp.behpath exp.name 'EMG_brusts_update'])%行为数据

load([exp.behpath exp.name 'EMG_brusts'])% 11-15行单独提取

load TL_ALL_include_EMG% in EMGfolder所有EMG数据

load([exp.behpath exp.name '_ALL_FuzzEn'])% fuzzy 数据
%% 1.1一个个验证RT更新是否有效， 999 为停止

% 获取所有显示器的位置
monitorPositions = get(0, 'MonitorPositions');

% 检查是否有多个屏幕
if size(monitorPositions, 1) < 2
    error('No second monitor detected');
end

% 获取第二个屏幕的位置
secondMonitor = monitorPositions(2, :);

% 创建figure窗口，并将其放置在第二个屏幕上
figure('Position', [secondMonitor(1), secondMonitor(2), 800, 600], 'WindowStyle', 'normal'); % 你可以根据需要调整宽度和高度


% 创建一个数组来记录每次用户输入
EMGonsite_Valid = 999*ones(6144, 1);

% 循环6144次
for i = 129:6144
    row = i;
    
    trial_to_plot = row; % 画EMG原始信号
    
    delay = 5; % 延迟60毫秒
    fs = 2000; % FDI 频率
    fs_bcp = 2000; % BCP 频率
    rt = AllBehaviour_EMGonsite(trial_to_plot, 6); % 反应时间
    FDIorBCP = AllBehaviour_EMGonsite(trial_to_plot, 4); % 1为FDI, 2为BCP
    evendt = AllBehaviour_EMGonsite(trial_to_plot, 7); % trial结束时间

    % 选择数据源并设置绘图参数
    if FDIorBCP == 1
        data_pre = abs([totEMG{trial_to_plot, 2}; totEMG{trial_to_plot, 3}(1:1000, :)]);
        ybound = [0, 2e-3];
        titlem = 'FDI - ';
    elseif FDIorBCP == 2
        data_pre = abs([totEMG_bcp{trial_to_plot, 2}; totEMG_bcp{trial_to_plot, 3}(1:1000, :)]);
        ybound = [0, 1e-3];
        titlem = 'BCP - ';
    end

    % 定义时间轴
    time_pre = (1:length(data_pre)) / fs;

    % 绘制EMG数据
    h = plot(time_pre, data_pre(:, 1), '-m');
    hold on;
    h2 = plot(time_pre, data_pre(:, 2), '-c');
    ylim(ybound);
    ylabel('Raw EMG');
    legend('left', 'right', 'left raw', 'right raw', 'AutoUpdate', 'off');

    % 添加时间线和标记
    xline(([0, 800, 1500, 2000, rt * 1000, evendt * 1000] + delay) / 1000, '--k', ...
          {'evidence on', 'minEvd0.8', 'DDL-1.5s', 'DDL-2s', 'RT', 'EVend'});

    % 绘制左手运动起点
    for j = 1:movement_counts_left(trial_to_plot)
        if movement_start_times_left{trial_to_plot}(j) < rt
            xline(movement_start_times_left{trial_to_plot}(j), '--m', 'LineWidth', 1.5);
        end
    end

    % 绘制右手运动起点
    for j = 1:movement_counts_right(trial_to_plot)
        if movement_start_times_right{trial_to_plot}(j) < rt
            xline(movement_start_times_right{trial_to_plot}(j), '--c', 'LineWidth', 1.5);
        end
    end

    % 设置标题
    title([titlem num2str(trial_to_plot)]);

    % 更新绘图
    drawnow;
    
    hold off;

    % 记录用户输入（1为有效，0为无效）
    prompt = 'Enter 1 for valid or 0 for invalid: ';
    userInput = input(prompt);
    if userInput==999
        break
    elseif isempty(userInput)
        EMGonsite_Valid(i) =1;

    else
    EMGonsite_Valid(i) = userInput;  % 将用户输入保存到数组中
    end

    % 显示当前输入
    disp(['Trial ', num2str(i), ': Input = ', num2str(userInput)]);
end
