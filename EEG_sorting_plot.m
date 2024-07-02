%%  3. EEG sortintg and plot
% Author: Liang Tong
% Date: 2/7/2024 

%% Toolbox requirements: 

clc
clear
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis'));% current folder

%% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);

%% Run analysis
eeglab


%% combining all data set
addpath(genpath(exp.finalpath));
load TL_ALL_include_EMG.mat% 包含所有行为数据包括EMG
%append 4 coloum
SL_b=999*ones(6144,1);%证据前0.188秒
SL_bb=999*ones(6144,1);% cue前0.188秒
RL_b=999*ones(6144,1);% 整段平均
RL_bb=999*ones(6144,1);% response 前0.188秒

for sub = exp.sub_id(1:end)
    [exp.finalpath 'Rejected_b_cICAri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.mat']

end

%% Sorting

