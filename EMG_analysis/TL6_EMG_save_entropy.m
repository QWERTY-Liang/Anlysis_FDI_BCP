%%  6.1 EMG analysis
% Author: Liang Tong
% Date: 12/7/2024 

%to update list
%1.function all the code to make main script simplier
%% Toolbox requirements: 
clc
clear all
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis'));% current folder

%% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);


%% 1. Load rl_BB pre-cue-188
% col 1: subject num
% col 2: block num
% col 3: contrast (high or low contrast)
% col 4: muscle (1=FDI, 2=BCP)
% col 5: trial outcome (1=correct, 2=error, 4=no response, 3 too early, 6 slow, 5 wrong muscle)
% col 6: RT in sec
% col 7: evshowtime =[]   % 0.2sec additional time after response
% col 8: participant response, 1 = left, 2 = right 0= on response
% col 9: correct response, 1 = left, 2 = right
% col 10: first tilt used for SSVEP
%证据前0.188秒 
% cue前0.188秒 col:12
% 整段平均 col:13
% response 前0.188秒 col:14
addpath(genpath(exp.finalpath));

load TL_ALL_include_EMG_updated.mat% EMG missing BCP totEMG

load TL_ALL_include_EMG% in EMGfolder


%%
    fs = 2000; % 假设采样频率为1000 Hz，根据你的实际情况调整
    low_cutoff = 50; % 低截止频率
    high_cutoff = 150; % 高截止频率
    filter_order = 4; % 滤波器阶数
    
    % 设计FIR带通滤波器
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
%%%需要调参
window_size=250;%越大约平滑，但是时间准确度变低
step_size=5;%越大数据点越密集，方便后续处理，但计算时间变长
delay=60;
fs=2000; % FDI 频率为2000
fs_bcp=2000;% BCP 频率为1926

    %% 1174 1175 1176  EMG wrong
    
% % Preallocate a temporary structure to store the results
FuzzEn_temp = struct('FuzzEn1', cell(6144, 1), 'FuzzEn2', cell(6144, 1), 'time_vector', cell(6144, 1));
N = 6144;
%parfor_progress(N); % Initialize 进度条
%Jeremy (2024). Progress monitor (progress bar) that works with parfor (https://www.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor-progress-bar-that-works-with-parfor), MATLAB Central File Exchange. Retrieved August 15, 2024.

parfor (trial_to_plot = 1:N,10)%最大核心数6
    % Initialize variables
    data_pre = [];
    filtered_data = [];
    filtered_data2 = [];
    time_vector=[];

    % Get the condition for FDIorBCP
    FDIorBCP = AllBehaviour_new(trial_to_plot, 4); % 1 for FDI, 2 for BCP

    % Prepare data based on the condition
    if FDIorBCP == 1
        data_pre = [totEMG{trial_to_plot, 2}; totEMG{trial_to_plot, 3}(1:1000, :)];
    elseif FDIorBCP == 2
        data_pre = [totEMG_bcp{trial_to_plot, 2}; totEMG_bcp{trial_to_plot, 3}(1:1000, :)];
    end

    % Check if data_pre is non-empty to avoid processing empty data
    if isempty(data_pre)|| any(isnan(data_pre), 'all') ||any(isinf(data_pre), 'all')
        disp(['The trial ' num2str(trial_to_plot) ' contains NaN or infinite values.']);
        continue;
    end

    time_pre = (1:length(data_pre)) / fs;

    % Separate data channels
    data = data_pre(:, 1); 
    data2 = data_pre(:, 2);

 % Apply filter notch
    data = filtfilt(notchFilter_50Hz, data);
    data2 = filtfilt(notchFilter_50Hz, data2);
    % Apply filter
    filtered_data = filtfilt(bp_filter, data);
    filtered_data2 = filtfilt(bp_filter, data2);

    % Compute Fuzzy Entropy using sliding window
    FuzzEn1 = TL6_sliding_window(filtered_data, window_size, step_size); % m=2; r=0.25
    FuzzEn2 = TL6_sliding_window(filtered_data2, window_size, step_size); % m=2; r=0.25

    % Time vector
    time_vector = ((1:length(FuzzEn1)) / length(FuzzEn1)) * (length(filtered_data) / fs);

    % Store results in FuzzEn_temp structure
    FuzzEn_temp(trial_to_plot).FuzzEn1 = FuzzEn1;
    FuzzEn_temp(trial_to_plot).FuzzEn2 = FuzzEn2;
    FuzzEn_temp(trial_to_plot).time_vector = time_vector;


    %parfor_progress; % Count计数进度条
end
%parfor_progress(0); % Clean up进度条
% Transfer results from FuzzEn_temp to FuzzEn


FuzzEn = cell(6144, 3);
for trial_to_plot = 1:6144
    FuzzEn{trial_to_plot, 1} = FuzzEn_temp(trial_to_plot).FuzzEn1;
    FuzzEn{trial_to_plot, 2} = FuzzEn_temp(trial_to_plot).FuzzEn2;
    FuzzEn{trial_to_plot, 3} = FuzzEn_temp(trial_to_plot).time_vector;
end

% Save as .mat file
save([exp.behpath exp.name '_ALL_FuzzEn'], 'FuzzEn','window_size','step_size');


% Example usage:
% 
% N = 100;
% parfor_progress(N); % Initialize
% parfor i=1:N
% pause(rand); % Replace with real code
% parfor_progress; % Count
% end
% parfor_progress(0); % Clean up

%% plot for testing
load([exp.behpath exp.name '_ALL_FuzzEn'])
% Selected rows to plot (e.g., rows 1, 2, 3, and 4)
selected_rows = 3000+[1, 2, 3, 4, 5, 6]; 

% Create a new figure for the subplots
figure;


% Loop through each selected row and plot
for i = 1:length(selected_rows)
    row = selected_rows(i);
    %
    trial_to_plot= row; %画EMG原始信号
    
    delay=0; %由于延迟所有EVon等推迟60毫秒，而真实rt要-60，EMG时表现为所有marker后移60，rt不变
fs=2000; % FDI 频率为2000
fs_bcp=2000;% BCP 频率为1926
rt=AllBehaviour_new(trial_to_plot,6);% matlab 记录的反应时间
FDIorBCP=AllBehaviour_new(trial_to_plot,4); %1为用了FDI, 2 为BCP
evendt=AllBehaviour_new(trial_to_plot,7);%trial结束时间
if FDIorBCP==1
data_pre=[totEMG{trial_to_plot, 2};totEMG{trial_to_plot, 3}(1:1000,:)];
ybound=[-0.1e-3,2e-3];
titlem=['FDI - '];
elseif FDIorBCP==2
data_pre=[totEMG_bcp{trial_to_plot, 2};totEMG_bcp{trial_to_plot, 3}(1:1000,:)];% 读取数据
ybound=[-0.05e-3,0.25e-3];
titlem=['BCP - '];
end
time_pre=(1:length(data_pre))/fs;
%

    % Extract data from FuzzEn
    data1 = FuzzEn{row, 1};
    data2 = FuzzEn{row, 2};
    time_vector = FuzzEn{row, 3};
    
    % Create a subplot for each selected row
    subplot(3, 2, i); % Arrange subplots in a 2x2 grid
    
    % Plot first column as red
    plot(time_vector, data1, 'r');
    hold on;
    % Plot second column as blue
    plot(time_vector, data2, 'b');
% Add labels and legend
    xlabel('Time (s)');
    ylabel('Fuzzy Entropy');
    

 % Plot additional data as green on right y-axis
    yyaxis right
    plot(time_pre, data_pre(:,1), '-m');
    plot(time_pre, data_pre(:,2), '-c');
    ylim(ybound)
    ylabel('Raw EMG');
legend('left', 'right','left raw','right raw','AutoUpdate', 'off');
xline(([0, 800, 1500,2000,rt*1000,evendt*1000]+delay)/1000, '--k', { 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s','RT','EVend'});
    
    subtitle([titlem num2str(trial_to_plot)]);
    
    hold off;
end
