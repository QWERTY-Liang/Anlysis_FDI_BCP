%%--------------------------------------------------------------------------------------------
%%  2. EEG analysis script
% Author: Liang Tong
% Date: 4/6/2024

%% Toolbox requirements: 
%Biosig: for reading BDF file

clc
clear
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis'));% current folder

%% Set experimental analysis parameters
%exp.sub_id = [1,2,3,4,5,6];
exp.sub_id = [1];
[exp] = TLBF1_setup(exp);

%% Run analysis
eeglab

%% STEP 0: Convert data to EEGlab format
% Convert data from bdf to .set and .fdt
for sub = exp.sub_id(1:end)

    TLBF2_convertData(sub, exp)

end
readlocs