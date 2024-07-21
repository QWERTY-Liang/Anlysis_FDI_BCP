%%  5. EEG sortintg and plot for pre evidence SL 
% Author: Liang Tong
% Date: 21/7/2024 

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
load TL_ALL_include_EMG.mat% 包含所有行为数据包括EMG
%append 4 coloum
SL_b=zeros(6144,1);%证据前0.188秒 col:11
SL_bb=zeros(6144,1);% cue前0.188秒 col:12
RL_b=zeros(6144,1);% 整段平均 col:13
RL_bb=zeros(6144,1);% response 前0.188秒 col:14
col11=0;
col12=0;%check number
col13=0;
col14=0;

for sub = exp.sub_id(1:end)
    filename=['Rejected_cICAriSL_' exp.filterLab 'aac' exp.name num2str(sub) '.mat'];
    load(filename)
    SL_b((sub-1)*1024 + 1:1024*sub)= bTrial_ind';
    col11=col11+length(bTrial_num);
    % 
    % filename=['Rejected_bb_cICAriSL_' exp.filterLab 'aac' exp.name num2str(sub) '.mat'];
    % load(filename)
    % SL_bb((sub-1)*1024 + 1:1024*sub)= bTrial_ind';
    % col12=col12+length(bTrial_num);

    filename=['Rejected_cICAriRL_' exp.filterLab 'aac' exp.name num2str(sub) '.mat'];
    load(filename)
    RL_b((sub-1)*1024 + 1:1024*sub)= bTrial_ind';
    col13=col13+length(bTrial_num);

    % filename=['Rejected_bb_cICAriRL_' exp.filterLab 'aac' exp.name num2str(sub) '.mat'];
    % load(filename)
    % RL_bb((sub-1)*1024 + 1:1024*sub)= bTrial_ind';
    % col14=col14+length(bTrial_num);

   % disp('finished')
end

if sum(SL_b)~=col11 || sum(SL_bb)~=col12 || sum(RL_b)~=col13 || sum(RL_bb)~=col14
    disp('ERROR: index not match!!!!')
else
    AllBehaviour_new=[AllBehaviour SL_b SL_bb  RL_b  RL_bb];
    % save as .mat file   
    save([exp.finalpath exp.name '_ALL_nobaselineEEG' ],'AllBehaviour_new','totEMG');


% Create variable names
    Varable_name={'subID' 'blocknum' 'contrast' 'muscle_used' 'perf' 'rt' 'evshowtime' 'respLR' 'corrLR' 'TiltSSVEP','SL_b','SL_bb','RL_b','RL_bb'};
    
   % Create tables with variable names
T_AllBehaviour_new = array2table(AllBehaviour_new, 'VariableNames', Varable_name);
% Write to a CSV file
writetable(T_AllBehaviour_new, [exp.finalpath exp.name '_ALL_ForR_nobaseline.csv']);
    disp(['Saved for updated' exp.name '_ALL_ for_R/ EMG' '.csv and .mat'])
end

%% 2. combine EEG

% case 1: combine EEG SL pre evidence (SL_bTLalldata.set)
% Initialize an empty ALLEEG structure
ALLEEG = [];
for sub = exp.sub_id(1:end)
     % % Load converted data
    EEG = pop_loadset([exp.finalpath 'csd_acICAri'  'SL_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, sub); % Store the dataset in ALLEEG
    disp(['File  ' num2str(sub) '  loaded successfully.']);
end
% Merge datasets
EEG = pop_mergeset(ALLEEG, 1:exp.sub_id(end), 0);

EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset

filename = ['nobaseline_SL_b' exp.name 'alldata'];
EEG = pop_saveset(EEG, filename, exp.finalpath);
% Display success message
disp(['File  ' filename '  saved successfully.']);

% sort matrix for case 1
AllBehaviour_SL_b=AllBehaviour_new;
for i=flip(1:length(AllBehaviour_new)) %倒序保证删除准确
    if AllBehaviour_new(i,11)==1
        AllBehaviour_SL_b(i,:)=[];
    end
end
 % save as .mat file   
    save([exp.finalpath exp.name '_AllBehaviour_SL_nobaseline' ],'AllBehaviour_SL_b');



    %% 2. combine EEG

% case 1: combine EEG SL pre evidence (SL_bTLalldata.set)
% Initialize an empty ALLEEG structure
ALLEEG = [];
for sub = exp.sub_id(1:end)
     % % Load converted data
    EEG = pop_loadset([exp.finalpath 'csd_acICAri'  'RL_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, sub); % Store the dataset in ALLEEG
    disp(['File  ' num2str(sub) '  loaded successfully.']);
end
% Merge datasets
EEG = pop_mergeset(ALLEEG, 1:exp.sub_id(end), 0);

EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset

filename = ['nobaseline_RL_b' exp.name 'alldata'];
EEG = pop_saveset(EEG, filename, exp.finalpath);
% Display success message
disp(['File  ' filename '  saved successfully.']);


load TL_ALL_nobaselineEEG.mat
% sort matrix for case 1
AllBehaviour_RL_bb= AllBehaviour_new;
for i=flip(1:length(AllBehaviour_new)) %倒序保证删除准确
    if AllBehaviour_new(i,13)==1
        AllBehaviour_RL_bb(i,:)=[];
    end
end
 % save as .mat file   
    save([exp.finalpath exp.name '_AllBehaviour_RL_nobaseline' ],'AllBehaviour_RL_bb');