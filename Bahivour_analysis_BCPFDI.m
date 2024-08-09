%%--------------------------------------------------------------------------------------------
%%  1. Bahivour analysis script
% Author: Liang Tong
% Date: 19/6/2024
%%
clc
clear
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis'));% current folder

%% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);

%% append all behaviour data for R/ feature analysis
% requiured output
% one .csv file for R analysis
% one .mt file include totEMG

% the behavioural data matrix should contain one row for each trial
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
%%%%%% col xxx: total EMG trail not ready


TLB1_AllBehaviour(exp);
% one .csv file for R analysis
% one .mt file include totEMG
% at EMG raw path
%%
% then modeling DDM, change the drift rate/boundary
addpath(genpath(exp.finalpath));

load TL_ALL_include_EMG_updated.mat% EMG missing BCP totEMG
