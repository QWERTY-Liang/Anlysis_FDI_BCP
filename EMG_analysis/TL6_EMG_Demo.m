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
% col 5: trial outcome (1=correct, 2=error, 4=no response, 3 too early, 6
% slow, 5 wrong muscle)  345 is useless data
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

%% plot sample trial
%partial brust: 3

%trial_to_plot=trial_to_plot+1;
trial_to_plot=7;
delay=60;
fs=2000; % FDI 频率为2000
fs_bcp=2000;% BCP 频率为1926
rt=AllBehaviour_new(trial_to_plot,6);% matlab 记录的反应时间
FDIorBCP=AllBehaviour_new(trial_to_plot,4); %1为用了FDI, 2 为BCP
evendt=AllBehaviour_new(trial_to_plot,7);%trial结束时间
% length(totEMG_bcp{trial_to_plot, 2})
% length(totEMG{trial_to_plot, 2})

figure;
time=(1:length(totEMG{trial_to_plot, 2})+length(totEMG{trial_to_plot, 3}))/fs;
time_bcp=(1:length(totEMG_bcp{trial_to_plot, 2})+length(totEMG_bcp{trial_to_plot, 3}))/fs_bcp;
% Plot totEMG
subplot(2, 1, 1);
plot(time, [totEMG{trial_to_plot, 2};totEMG{trial_to_plot, 3}]);
ylim([-2e-3,2e-3])
%xlim([0,2.5])
hold on;
xline(([0, 800, 1500,2000,rt*1000,evendt*1000]+delay)/1000, '--r', { 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s','RT','EVend'});
subtitle(['FDI - ' num2str(trial_to_plot)]);
hold off

% Plot totEMG_bcp
subplot(2, 1, 2);
plot(time_bcp, [totEMG_bcp{trial_to_plot, 2};totEMG_bcp{trial_to_plot, 3}]);
subtitle(['BCP - ' num2str(trial_to_plot)]);
%xlim([0,2.5])
ylim([-2e-3,2e-3])
xline(([0, 800, 1500,2000,rt*1000,evendt*1000]+delay)/1000, '--r', { 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s','RT','EVend'});



%% Fuzzy entropy test

 % 设计带通滤波器
    fs = 2000; % 假设采样频率为1000 Hz，根据你的实际情况调整
    low_cutoff = 20; % 低截止频率
    high_cutoff = 250; % 高截止频率
    filter_order = 4; % 滤波器阶数
    
    % 设计FIR带通滤波器
    bp_filter = designfilt('bandpassfir', 'FilterOrder', filter_order, ...
                           'CutoffFrequency1', low_cutoff, 'CutoffFrequency2', high_cutoff, ...
                           'SampleRate', fs);

window_size=150;%越大约平滑，但是时间准确度变低
step_size=1;%越大数据点越密集，方便后续处理，但计算时间变长

if FDIorBCP==1
data_pre=[totEMG{trial_to_plot, 2};totEMG{trial_to_plot, 3}(1:1000,:)];

titlem=['FDI - '];
elseif FDIorBCP==2
data_pre=[totEMG_bcp{trial_to_plot, 2};totEMG_bcp{trial_to_plot, 3}(1:1000,:)];% 读取数据
titlem=['BCP - '];
end
time_pre=(1:length(data_pre))/fs;



data=data_pre(:,1); % delete offest
data2=data_pre(:,2);
%滤波
filtered_data = filtfilt(bp_filter, data);
filtered_data2 = filtfilt(bp_filter, data2);
%plot(filtered_data )
output = TL6_sliding_window(filtered_data, window_size, step_size);%m=2;r=0.25;%应用滑动渐变窗，计算模糊熵
output2 = TL6_sliding_window(filtered_data2, window_size, step_size);%m=2;r=0.25;%应用滑动渐变窗，计算模糊熵
t_output =(   (1:length(output))/length(output)  ) *(length(filtered_data)/fs);




figure
subplot(2, 1, 1);
plot(t_output,output, t_output,output2)
xlim([0,2.5])
hold on;
xline(([0, 800, 1500,2000,rt*1000,evendt*1000]+delay)/1000, '--r', { 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s','RT','EVend'});
subtitle([titlem num2str(trial_to_plot)]);
hold off
subplot(2, 1, 2);
plot(time_pre, data_pre);
ylim([-2e-3,2e-3])
xlim([0,2.5])
hold on;
xline(([0, 800, 1500,2000,rt*1000,evendt*1000]+delay)/1000, '--r', { 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s','RT','EVend'});
subtitle([titlem num2str(trial_to_plot)]);
hold off

