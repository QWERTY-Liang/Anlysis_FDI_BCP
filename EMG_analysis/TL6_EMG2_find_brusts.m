%%  6.1 EMG analysis
% Author: Liang Tong
% Date: 12/7/2024 

%to update list
%1.function all the code to make main script simplier
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

%load TL_ALL_include_EMG_updated.mat% EMG missing BCP totEMG

load TL_ALL_include_EMG_1% in EMGfolder

load([exp.behpath exp.name '_ALL_FuzzEn_processed']) % sep 5 2024 version
%%
% Parameters
smoothing_window = 5; % Window size for smoothing (adjust as needed)
% threshold_factor_L = 4; % Factor to determine the threshold based on the mean of the data
% threshold_factor_R = 1.7;% 左边信号质量都差，所以右边可以更严格阈值
baselinewindow= 100;%调节平均熵阈值的长度

% Initialize arrays to store movement information for both hands
movement_counts_left = zeros(size(FuzzEn, 1), 1); % Number of movements detected for the left hand
movement_counts_right = zeros(size(FuzzEn, 1), 1); % Number of movements detected for the right hand
movement_start_times_left = cell(size(FuzzEn, 1), 1); % Cell array to store start times for the left hand
%movement_end_times_left = cell(size(FuzzEn, 1), 1); % Cell array to store end times for the left hand
movement_start_times_right = cell(size(FuzzEn, 1), 1); % Cell array to store start times for the right hand
%movement_end_times_right = cell(size(FuzzEn, 1), 1); % Cell array to store end times for the right hand

for trial_to_plot = 1:size(FuzzEn, 1)

%%%两种肌肉不同阈值
FDIorBCP=AllBehaviour(trial_to_plot,4); %1为用了FDI, 2 为BCP
if FDIorBCP==4
threshold_factor_L = 1.3; % Factor to determine the threshold based on the mean of the data
threshold_factor_R = 1.4;% 左边信号质量都差，所以右边可以更严格阈值

elseif FDIorBCP==44
threshold_factor_L = 1.2; % Factor to determine the threshold based on the mean of the data
threshold_factor_R = 1.3;% 左边信号质量都差，所以右边可以更严格阈值

end

    % Extract the entropy sequences for both hands
    entropy_seq_left = FuzzEn{trial_to_plot, 1};
    entropy_seq_right = FuzzEn{trial_to_plot, 2};
        % Check if data_pre is non-empty to avoid processing empty data
    if isempty(entropy_seq_left)|| any(isnan(entropy_seq_left), 'all') ||any(isinf(entropy_seq_left), 'all')||isempty(entropy_seq_right)|| any(isnan(entropy_seq_right), 'all') ||any(isinf(entropy_seq_right), 'all')
        disp(['The trial ' num2str(trial_to_plot) ' contains NaN or infinite values.']);
        continue;
    end
    % Smooth the entropy sequences to reduce noise
    % smoothed_entropy_left = smooth(entropy_seq_left, smoothing_window);
    % smoothed_entropy_right = smooth(entropy_seq_right, smoothing_window);
    %method 2 -mean then abs
    smoothed_entropy_left = abs(-entropy_seq_left-mean(entropy_seq_left(1:baselinewindow)));
    smoothed_entropy_right = abs(entropy_seq_right-mean(entropy_seq_right(1:baselinewindow)));
    
    % Determine thresholds based on the mean of the smoothed entropy
    threshold_left = mean(smoothed_entropy_left) * threshold_factor_L;% method 1 is mean
    threshold_right = mean(smoothed_entropy_right) * threshold_factor_R;% method 1 is mean
    if threshold_right<0.16%这里只是处理前信号使用
        threshold_right=0.16;%如果信号很稳定，设置最低阈值

    end
        if threshold_left<0.15
       threshold_left=0.15;
    end


    % Find indices where the smoothed entropy crosses the threshold
    above_threshold_left = smoothed_entropy_left > threshold_left;
    above_threshold_right = smoothed_entropy_right > threshold_right;
    
    % Identify start and end points of all segments above the threshold for the left hand
    starts_left = find(diff([0; above_threshold_left]) == 1); % Points where it starts to be above threshold
    ends_left = find(diff([above_threshold_left; 0]) == -1); % Points where it drops back below threshold
    
    % Identify start and end points of all segments above the threshold for the right hand
    starts_right = find(diff([0; above_threshold_right]) == 1); % Points where it starts to be above threshold
    ends_right = find(diff([above_threshold_right; 0]) == -1); % Points where it drops back below threshold
    
    % Record the number of movements detected for each hand
    movement_counts_left(trial_to_plot) = length(starts_left);
    movement_counts_right(trial_to_plot) = length(starts_right);
    
    % Record the start and end times for each detected segment
    movement_start_times_left{trial_to_plot} = FuzzEn{trial_to_plot, 3}(starts_left);
    %movement_end_times_left{trial_to_plot} = FuzzEn{trial_to_plot, 3}(ends_left);
    movement_start_times_right{trial_to_plot} = FuzzEn{trial_to_plot, 3}(starts_right);
    %movement_end_times_right{trial_to_plot} = FuzzEn{trial_to_plot, 3}(ends_right);
end
% 检查可用trial数
[sum((movement_counts_right==1 | movement_counts_right==2));sum(movement_counts_left==1 | movement_counts_left==2);sum(movement_counts_right==0);sum(movement_counts_left==0)]
% Save as .mat file
save([exp.behpath exp.name 'EMG_brusts'], 'movement_counts_left','movement_counts_right','movement_start_times_left','movement_start_times_right');

% % Display the results for each trial
% for trial_to_plot = 1:size(FuzzEn, 1)
%     fprintf('Trial %d:\n', trial_to_plot);
%     fprintf('Number of Movements Detected (Left Hand): %d\n', movement_counts_left(trial_to_plot));
%     if movement_counts_left(trial_to_plot) > 0
%         fprintf('Start Times (Left Hand): %s\n', num2str(movement_start_times_left{trial_to_plot}'));
%         fprintf('End Times (Left Hand): %s\n', num2str(movement_end_times_left{trial_to_plot}'));
%     end
%     fprintf('Number of Movements Detected (Right Hand): %d\n', movement_counts_right(trial_to_plot));
%     if movement_counts_right(trial_to_plot) > 0
%         fprintf('Start Times (Right Hand): %s\n', num2str(movement_start_times_right{trial_to_plot}'));
%         fprintf('End Times (Right Hand): %s\n', num2str(movement_end_times_right{trial_to_plot}'));
% %     end
% %     fprintf('\n');
% % end
% 单独作画


% % Optional: plot the results for a specific trial
% trial_to_plot = 4; % Choose a specific trial to plot
% rt=AllBehaviour_new(trial_to_plot,6);% matlab 记录的反应时间
% %读取源数据加命名
% if FDIorBCP==1
% data_pre=[totEMG{trial_to_plot, 2};totEMG{trial_to_plot, 3}(1:1000,:)];
% 
% titlem=['FDI - '];
% elseif FDIorBCP==2
% data_pre=[totEMG_bcp{trial_to_plot, 2};totEMG_bcp{trial_to_plot, 3}(1:1000,:)];% 读取数据
% titlem=['BCP - '];
% end
% time_pre=(1:length(data_pre))/fs;
% 
% figure;
% % Plot left-hand entropy with blue color
% %plot(FuzzEn{trial_to_plot, 3}, FuzzEn{trial_to_plot, 1}, 'b'); % smooth
% plot(FuzzEn{trial_to_plot, 3}, abs(FuzzEn{trial_to_plot, 1}-mean(FuzzEn{trial_to_plot, 1}(1:100))), 'b'); %method 2 mean abs
% hold on;
% % Plot right-hand entropy with red color
% plot(FuzzEn{trial_to_plot, 3}, abs(FuzzEn{trial_to_plot, 2}-mean(FuzzEn{trial_to_plot, 2}(1:100))), 'r'); 
% 
% % Plot detected start and end points for the left hand with green and magenta lines
% for i = 1:movement_counts_left(trial_to_plot)
%     xline(movement_start_times_left{trial_to_plot}(i), 'g'); % Start point
%     xline(movement_end_times_left{trial_to_plot}(i), 'm'); % End point
% end
% 
% % Plot detected start and end points for the right hand with cyan and yellow lines
% for i = 1:movement_counts_right(trial_to_plot)
%     xline(movement_start_times_right{trial_to_plot}(i), 'c'); % Start point
%     xline(movement_end_times_right{trial_to_plot}(i), 'y'); % End point
% end
% xline(rt,'k', 'LineWidth', 2)
% xlabel('Time (s)');
% ylabel('Fuzzy Entropy');
% legend('Left Hand Entropy', 'Right Hand Entropy', 'Left Start Points', 'Left End Points', 'Right Start Points', 'Right End Points', 'RT');
% title(['FuzzEn for Trial ' num2str(trial_to_plot)]);
% hold off;


%% 多幅图同时画

selected_rows =10+[1, 2, 3, 4, 5, 6,7,8]; 

% Create a new figure for the subplots
figure;
% Loop through each selected row and plot
for i = 1:length(selected_rows)
    row = selected_rows(i);
    %
    trial_to_plot= row; %画EMG原始信号
    
    delay=5; %由于延迟所有EVon等推迟60毫秒，而真实rt要-60，EMG时表现为所有marker后移60，rt不变
fs=2000; % FDI 频率为2000
fs_bcp=2000;% BCP 频率为1926
rt=AllBehaviour(trial_to_plot,6);% matlab 记录的反应时间
FDIorBCP=AllBehaviour(trial_to_plot,4); %4为用了FDI, 44 为BCP
evendt=AllBehaviour(trial_to_plot,7);%trial结束时间

if FDIorBCP==4
data_pre=abs([totEMG{trial_to_plot, 2};totEMG{trial_to_plot, 3}(1:1000,:)]);
ybound=[0,2e-3];
titlem=['FDI - '];
elseif FDIorBCP==44
data_pre=abs([totEMG_bcp{trial_to_plot, 2};totEMG_bcp{trial_to_plot, 3}(1:1000,:)]);% 读取数据
ybound=[0,1e-3];
titlem=['BCP - '];
end
time_pre=(1:length(data_pre))/fs;
%

    % Extract data from FuzzEn
    data1 = abs((FuzzEn{row, 1}-mean(FuzzEn{row, 1}(1:baselinewindow))));
    data2 = abs((FuzzEn{row, 2}-mean(FuzzEn{row, 2}(1:baselinewindow))));
    time_vector = FuzzEn{row, 3};
    
    % Create a subplot for each selected row
    subplot(4, 2, i); % Arrange subplots in a 2x2 grid
    
    % Plot first column as red
    plot(time_vector, data1, 'r', 'LineWidth', 0.5);
    hold on;
    % Plot second column as blue
    plot(time_vector, data2, 'b', 'LineWidth', 0.5);

    % Plot detected start and end points for the left hand with green and magenta lines

%xline(rt,'k', 'LineWidth', 2)
xlabel('Time (s)');
ylabel('Fuzzy Entropy');
legend('Left Hand Entropy', 'Right Hand Entropy')

 % Plot additional data as green on right y-axis
    yyaxis right
    h=plot(time_pre, data_pre(:,1), '-m');
    h2=plot(time_pre, data_pre(:,2), '-c');
    % alpha(h,.5)
    % alpha(h2,.5)
    ylim(ybound)
    ylabel('Raw EMG');
legend('left', 'right','left raw','right raw','AutoUpdate', 'off');
xline(([0, 800, 1500,2000,rt*1000,evendt*1000]+delay)/1000, '--k', { 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s','RT','EVend'});
   for i = 1:movement_counts_left(trial_to_plot)
    xline(movement_start_times_left{trial_to_plot}(i), '--m', 'LineWidth', 1.5); % Start point
    %xline(movement_end_times_left{trial_to_plot}(i), '--r', 'LineWidth', 1.2); % End point
end

% Plot detected start and end points for the right hand with cyan and yellow lines
for i = 1:movement_counts_right(trial_to_plot)
    xline(movement_start_times_right{trial_to_plot}(i), '--c', 'LineWidth', 1.5); % Start point
    %xline(movement_end_times_right{trial_to_plot}(i), '--b', 'LineWidth', 1.2); % End point
end 
    subtitle([titlem num2str(trial_to_plot)]);
    
    hold off;
end

%% 检查可用trial数
sum((movement_counts_right==1 | movement_counts_right==2))
sum(movement_counts_left==1 | movement_counts_left==2)


sum(movement_counts_right==0)
sum(movement_counts_left==0)


