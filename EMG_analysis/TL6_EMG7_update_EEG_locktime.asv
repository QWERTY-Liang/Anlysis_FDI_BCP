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
%col 15 是是否更新
%col 16 手检后EMG_onsite 时间
addpath(genpath(exp.finalpath));

% load TL_ALL_include_EMG_updated.mat% EMG missing BCP totEMG
% 
% load TL_ALL_include_EMG% in EMGfolder

load([exp.behpath exp.name 'EMG_brusts_update'])

load([exp.behpath exp.name 'EMG_Vaild_check'])%trial目检有效性
load([exp.behpath exp.name 'EMGonsite_Lock_time'])%更新后时间

%%

