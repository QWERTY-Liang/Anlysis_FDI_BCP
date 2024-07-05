%%  4.1 EEG sortintg and plot for pre cue SL  LRP 
% Author: Liang Tong
% Date: 4/7/2024 

%to update list
%1.function all the code to make main script simplier
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
%% 1. combining all behaviour data set

addpath(genpath(exp.finalpath));