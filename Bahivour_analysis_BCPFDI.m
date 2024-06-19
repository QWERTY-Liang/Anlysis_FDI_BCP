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
%exp.sub_id = [1,2,3,4,5,6];
exp.sub_id = [1];
[exp] = TLBF1_setup(exp);

%% accuracy analysis
%1. indivial 
%2. overall