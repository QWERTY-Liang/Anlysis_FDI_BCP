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

load([exp.behpath exp.name '_ALL_FuzzEn'])
%%
% Parameters
smoothing_window = 5; % Window size for smoothing (adjust as needed)
threshold_factor_L = 1.2; % Factor to determine the threshold based on the mean of the data
threshold_factor_R = 1.5;% 左边信号质量都差，所以右边可以更严格阈值

% Initialize arrays to store movement information for both hands
movement_counts_left = zeros(size(FuzzEn, 1), 1); % Number of movements detected for the left hand
movement_counts_right = zeros(size(FuzzEn, 1), 1); % Number of movements detected for the right hand
movement_start_times_left = cell(size(FuzzEn, 1), 1); % Cell array to store start times for the left hand
movement_end_times_left = cell(size(FuzzEn, 1), 1); % Cell array to store end times for the left hand
movement_start_times_right = cell(size(FuzzEn, 1), 1); % Cell array to store start times for the right hand
movement_end_times_right = cell(size(FuzzEn, 1), 1); % Cell array to store end times for the right hand

for trial_to_plot = 1:size(FuzzEn, 1)
    % Extract the entropy sequences for both hands
    entropy_seq_left = FuzzEn{trial_to_plot, 1};
    entropy_seq_right = FuzzEn{trial_to_plot, 2};
    
    % Smooth the entropy sequences to reduce noise
    smoothed_entropy_left = smooth(entropy_seq_left, smoothing_window);
    smoothed_entropy_right = smooth(entropy_seq_right, smoothing_window);
    
    % Determine thresholds based on the mean of the smoothed entropy
    threshold_left = mean(smoothed_entropy_left) * threshold_factor_L;
    threshold_right = mean(smoothed_entropy_right) * threshold_factor_R;
    
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
    movement_end_times_left{trial_to_plot} = FuzzEn{trial_to_plot, 3}(ends_left);
    movement_start_times_right{trial_to_plot} = FuzzEn{trial_to_plot, 3}(starts_right);
    movement_end_times_right{trial_to_plot} = FuzzEn{trial_to_plot, 3}(ends_right);
end

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
%     end
%     fprintf('\n');
% end

% Optional: plot the results for a specific trial
trial_to_plot = 6000; % Choose a specific trial to plot
rt=AllBehaviour_new(trial_to_plot,6);% matlab 记录的反应时间
figure;
% Plot left-hand entropy with blue color
plot(FuzzEn{trial_to_plot, 3}, FuzzEn{trial_to_plot, 1}, 'b'); 
hold on;
% Plot right-hand entropy with red color
plot(FuzzEn{trial_to_plot, 3}, FuzzEn{trial_to_plot, 2}, 'r'); 

% Plot detected start and end points for the left hand with green and magenta lines
for i = 1:movement_counts_left(trial_to_plot)
    xline(movement_start_times_left{trial_to_plot}(i), 'g'); % Start point
    xline(movement_end_times_left{trial_to_plot}(i), 'm'); % End point
end

% Plot detected start and end points for the right hand with cyan and yellow lines
for i = 1:movement_counts_right(trial_to_plot)
    xline(movement_start_times_right{trial_to_plot}(i), 'c'); % Start point
    xline(movement_end_times_right{trial_to_plot}(i), 'y'); % End point
end
xline(rt,'k', 'LineWidth', 2)
xlabel('Time (s)');
ylabel('Fuzzy Entropy');
legend('Left Hand Entropy', 'Right Hand Entropy', 'Left Start Points', 'Left End Points', 'Right Start Points', 'Right End Points', 'RT');
title(['FuzzEn for Trial ' num2str(trial_to_plot)]);
hold off;
