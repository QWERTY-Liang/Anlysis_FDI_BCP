%%  3.2 EEG sortintg and plot for pre cue SL 
% Author: Liang Tong
% Date: 2/7/2024 

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
% load TL_ALL_include_EMG.mat% 包含所有行为数据包括EMG
% %append 4 coloum
% SL_b=zeros(6144,1);%证据前0.188秒 col:11
% SL_bb=zeros(6144,1);% cue前0.188秒 col:12
% RL_b=zeros(6144,1);% 整段平均 col:13
% RL_bb=zeros(6144,1);% response 前0.188秒 col:14
% col11=0;
% col12=0;%check number
% col13=0;
% col14=0;
% 
% for sub = exp.sub_id(1:end)
%     filename=['Rejected_b_cICAriSL_' exp.filterLab 'aac' exp.name num2str(sub) '.mat'];
%     load(filename)
%     SL_b((sub-1)*1024 + 1:1024*sub)= bTrial_ind';
%     col11=col11+length(bTrial_num);
%     filename=['Rejected_bb_cICAriSL_' exp.filterLab 'aac' exp.name num2str(sub) '.mat'];
%     load(filename)
%     SL_bb((sub-1)*1024 + 1:1024*sub)= bTrial_ind';
%     col12=col12+length(bTrial_num);
%     filename=['Rejected_b_cICAriRL_' exp.filterLab 'aac' exp.name num2str(sub) '.mat'];
%     load(filename)
%     RL_b((sub-1)*1024 + 1:1024*sub)= bTrial_ind';
%     col13=col13+length(bTrial_num);
%     filename=['Rejected_bb_cICAriRL_' exp.filterLab 'aac' exp.name num2str(sub) '.mat'];
%     load(filename)
%     RL_bb((sub-1)*1024 + 1:1024*sub)= bTrial_ind';
%     col14=col14+length(bTrial_num);
% 
%    % disp('finished')
% end
% 
% if sum(SL_b)~=col11 || sum(SL_bb)~=col12 || sum(RL_b)~=col13 || sum(RL_bb)~=col14
%     disp('ERROR: index not match!!!!')
% else
%     AllBehaviour_new=[AllBehaviour SL_b SL_bb  RL_b  RL_bb];
%     % save as .mat file   
%     save([exp.finalpath exp.name '_ALL_include_EMG_updated' ],'AllBehaviour_new','totEMG');
% 
% 
% % Create variable names
%     Varable_name={'subID' 'blocknum' 'contrast' 'muscle_used' 'perf' 'rt' 'evshowtime' 'respLR' 'corrLR' 'TiltSSVEP','SL_b','SL_bb','RL_b','RL_bb'};
% 
%    % Create tables with variable names
% T_AllBehaviour_new = array2table(AllBehaviour_new, 'VariableNames', Varable_name);
% % Write to a CSV file
% writetable(T_AllBehaviour_new, [exp.finalpath exp.name '_ALL_ForR_updated.csv']);
%     disp(['Saved for updated' exp.name '_ALL_ for_R/ EMG' '.csv and .mat'])
% end

%% 2. combine EEG

% case 1: combine EEG SL pre evidence (SL_bTLalldata.set)
% Initialize an empty ALLEEG structure
ALLEEG = [];
for sub = exp.sub_id(1:end)
     % % Load converted data
    EEG = pop_loadset([exp.finalpath 'csd_abb_cICAri'  'SL_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, sub); % Store the dataset in ALLEEG
    disp(['File  ' num2str(sub) '  loaded successfully.']);
end
% Merge datasets
EEG = pop_mergeset(ALLEEG, 1:exp.sub_id(end), 0);

EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset

filename = ['SL_bb' exp.name 'alldata'];
EEG = pop_saveset(EEG, filename, exp.finalpath);
% Display success message
disp(['File  ' filename '  saved successfully.']);


load TL_ALL_include_EMG_updated.mat
% sort matrix for case 1
AllBehaviour_SL_bb= AllBehaviour_new;
for i=flip(1:length(AllBehaviour_new)) %倒序保证删除准确
    if AllBehaviour_new(i,12)==1
        AllBehaviour_SL_bb(i,:)=[];
    end
end
 % save as .mat file   
    save([exp.finalpath exp.name '_AllBehaviour_SL_bb' ],'AllBehaviour_SL_bb');


%% 3. Sorting v%% pre cue SL
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
 EEG = pop_loadset([exp.finalpath 'SL_bbTLalldata.set']);
 load TL_AllBehaviour_SL_bb.mat 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.1 correct vs wrong
% Define the selected trials vector
channel_to_plot = 19;

selected_trials_1 = zeros(length(AllBehaviour_SL_bb),1);  % example vector
selected_trials_2 = zeros(length(AllBehaviour_SL_bb),1);  % example vector
count=0;
% switch
%     case 'CW_SL_b' %correct vs wrong; pre cue SL
        for i=1:length(AllBehaviour_SL_bb)
            %select trials
            if AllBehaviour_SL_bb(i,8)==AllBehaviour_SL_bb(i,9)
                selected_trials_1(i)=1;
            elseif AllBehaviour_SL_bb(i,8)~=AllBehaviour_SL_bb(i,9)
                selected_trials_2(i)=1;
            end
            %exclude invalide(too early or wrong muscle)
            if AllBehaviour_SL_bb(i,5)==3 || AllBehaviour_SL_bb(i,5)==5
                selected_trials_1(i)=0;%change back to unselected
                selected_trials_2(i)=0;
                count=count+1;
            end
        end
% end

% Extract selected trials---for correct
selected_idx_1 = find(selected_trials_1 == 1);
EEG_selected_2 = pop_select(EEG, 'trial', selected_idx_1);
% Compute the ERP
ERP_correct = mean(EEG_selected_2.data, 3);  % averaging across the third dimension (trials)
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_correct = EEG_selected_2.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% Extract selected trials---for wrong
selected_idx_wrong = find(selected_trials_2 == 1);
EEG_selected_wrong = pop_select(EEG, 'trial', selected_idx_wrong);
% Compute the ERP
ERP_wrong = mean(EEG_selected_wrong.data, 3);  % averaging across the third dimension (trials)
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_wrong = EEG_selected_wrong.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% plot ERPs
figure;
plot(time_vector_correct, ERP_correct(channel_to_plot, :));
hold on;
plot(time_vector_wrong, ERP_wrong(channel_to_plot, :));
legend('correct','wrong')
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(['ERP for Channel ', num2str(channel_to_plot)]);
hold off;
% Add vertical lines at specified time points
hold on;
xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s'});
hold off;

%% plot scalp plot ----- correct
num_plots = 10;
time_points = linspace(-1200, 1500, num_plots);  % in milliseconds

% Find the indices of the defined time points
time_indices = arrayfun(@(t) find(EEG_selected_2.times >= t, 1), time_points);

% Plot scalp topographies
figure;
for i = 1:num_plots
    subplot(2, 5, i);  % create a 2x4 subplot
    topoplot(ERP_correct(:, time_indices(i)), EEG_selected_2.chanlocs, 'maplimits', [-max(abs(ERP_correct(:))), max(abs(ERP_correct(:)))]);
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end

% plot scalp plot ----- wrong
num_plots = 10;
time_points = linspace(-1200, 1500, num_plots);  % in milliseconds

% Find the indices of the defined time points
time_indices = arrayfun(@(t) find(EEG_selected_wrong.times >= t, 1), time_points);

% Plot scalp topographies
figure;
for i = 1:num_plots
    subplot(2, 5, i);  % create a 2x4 subplot
    topoplot(ERP_wrong(:, time_indices(i)), EEG_selected_wrong.chanlocs, 'maplimits', [-max(abs(ERP_wrong(:))), max(abs(ERP_wrong(:)))]);
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3.2 high Vs low
% Define the selected trials vector
channel_to_plot = 19;

selected_trials_1 = zeros(length(AllBehaviour_SL_bb),1);  % example vector
selected_trials_2 = zeros(length(AllBehaviour_SL_bb),1);  % example vector
count=0;

for i=1:length(AllBehaviour_SL_bb)

    if AllBehaviour_SL_bb(i,3)==0.14
        selected_trials_1(i)=1;
    elseif AllBehaviour_SL_bb(i,3)==0.07
        selected_trials_2(i)=1;
    end


    %exclude invalide(too early or wrong muscle)
    if AllBehaviour_SL_bb(i,5)==3 || AllBehaviour_SL_bb(i,5)==5
        selected_trials_1(i)=0;%change back to unselected
        selected_trials_2(i)=0;
        count=count+1;
    end
end

% Extract selected trials---for correct
selected_idx_1 = find(selected_trials_1 == 1);
EEG_selected_2 = pop_select(EEG, 'trial', selected_idx_1);
% Compute the ERP
ERP_correct = mean(EEG_selected_2.data, 3);  % averaging across the third dimension (trials)
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_correct = EEG_selected_2.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% Extract selected trials---for wrong
selected_idx_wrong = find(selected_trials_2 == 1);
EEG_selected_wrong = pop_select(EEG, 'trial', selected_idx_wrong);
% Compute the ERP
ERP_wrong = mean(EEG_selected_wrong.data, 3);  % averaging across the third dimension (trials)
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_wrong = EEG_selected_wrong.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% plot ERPs
figure;
plot(time_vector_correct, ERP_correct(channel_to_plot, :));
hold on;
plot(time_vector_wrong, ERP_wrong(channel_to_plot, :));
legend('High(0.14)','Low(0.07)')
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(['ERP for Channel ', num2str(channel_to_plot)]);
hold off;
% Add vertical lines at specified time points
hold on;
xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s'});
hold off;

%% plot scalp plot ----- High
num_plots = 10;
time_points = linspace(-1200, 1500, num_plots);  % in milliseconds

% Find the indices of the defined time points
time_indices = arrayfun(@(t) find(EEG_selected_2.times >= t, 1), time_points);

% Plot scalp topographies
figure;
for i = 1:num_plots
    subplot(2, 5, i);  % create a 2x4 subplot
    topoplot(ERP_correct(:, time_indices(i)), EEG_selected_2.chanlocs, 'maplimits', [-max(abs(ERP_correct(:))), max(abs(ERP_correct(:)))]);
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end

% plot scalp plot ----- Low
num_plots = 10;
time_points = linspace(-1200, 1500, num_plots);  % in milliseconds

% Find the indices of the defined time points
time_indices = arrayfun(@(t) find(EEG_selected_wrong.times >= t, 1), time_points);
% Plot scalp topographies
figure;
for i = 1:num_plots
    subplot(2, 5, i);  % create a 2x4 subplot
    topoplot(ERP_wrong(:, time_indices(i)), EEG_selected_wrong.chanlocs, 'maplimits', [-max(abs(ERP_wrong(:))), max(abs(ERP_wrong(:)))]);
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.3 FDI vs BCP
% Define the selected trials vector
channel_to_plot = 19;

selected_trials_1 = zeros(length(AllBehaviour_SL_bb),1);  % example vector
selected_trials_2 = zeros(length(AllBehaviour_SL_bb),1);  % example vector
count=0;

for i=1:length(AllBehaviour_SL_bb)


    %FDI and BCP only
    % if AllBehaviour_SL_bb(i,4)==1
    %     selected_trials_1(i)=1;
    % elseif AllBehaviour_SL_bb(i,4)==2
    %     selected_trials_2(i)=1;
    % end
%%%%%%%%%%%%%%%%%%%%%change here
    %fdi correct and BCP correct (==) or wrong(~=); high=0.14 low 0.07
    if AllBehaviour_SL_bb(i,4)==1 &&AllBehaviour_SL_bb(i,8)~=AllBehaviour_SL_bb(i,9)&& AllBehaviour_SL_bb(i,3)==0.14
        selected_trials_1(i)=1;
    elseif AllBehaviour_SL_bb(i,4)==2 &&AllBehaviour_SL_bb(i,8)~=AllBehaviour_SL_bb(i,9)&& AllBehaviour_SL_bb(i,3)==0.14
        selected_trials_2(i)=1;
    end


    %exclude invalide(too early or wrong muscle)
    if AllBehaviour_SL_bb(i,5)==3 || AllBehaviour_SL_bb(i,5)==5
        selected_trials_1(i)=0;%change back to unselected
        selected_trials_2(i)=0;
        count=count+1;
    end
end

% Extract selected trials---for correct
selected_idx_1 = find(selected_trials_1 == 1);
EEG_selected_2 = pop_select(EEG, 'trial', selected_idx_1);
% Compute the ERP
ERP_correct = mean(EEG_selected_2.data, 3);  % averaging across the third dimension (trials)
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_correct = EEG_selected_2.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% Extract selected trials---for wrong
selected_idx_wrong = find(selected_trials_2 == 1);
EEG_selected_wrong = pop_select(EEG, 'trial', selected_idx_wrong);
% Compute the ERP
ERP_wrong = mean(EEG_selected_wrong.data, 3);  % averaging across the third dimension (trials)
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_wrong = EEG_selected_wrong.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% plot ERPs
figure;
plot(time_vector_correct, ERP_correct(channel_to_plot, :));
hold on;
plot(time_vector_wrong, ERP_wrong(channel_to_plot, :));
legend('FDI','BCP')
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(['ERP for Channel ', num2str(channel_to_plot)]);
hold off;
% Add vertical lines at specified time points
hold on;
xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s'});
hold off;

%% plot scalp plot ----- FDI
num_plots = 10;
time_points = linspace(-600, 750, num_plots);  % in milliseconds

% Find the indices of the defined time points
time_indices = arrayfun(@(t) find(EEG_selected_2.times >= t, 1), time_points);

% Plot scalp topographies
figure;
for i = 1:num_plots
    subplot(2, 5, i);  % create a 2x4 subplot
    topoplot(ERP_correct(:, time_indices(i)), EEG_selected_2.chanlocs, 'maplimits', [-max(abs(ERP_correct(:))), max(abs(ERP_correct(:)))]);
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end

% plot scalp plot ----- BCP
num_plots = 10;
time_points = linspace(-600, 750, num_plots);  % in milliseconds

% Find the indices of the defined time points
time_indices = arrayfun(@(t) find(EEG_selected_wrong.times >= t, 1), time_points);

% Plot scalp topographies
figure;
for i = 1:num_plots
    subplot(2, 5, i);  % create a 2x4 subplot
    topoplot(ERP_wrong(:, time_indices(i)), EEG_selected_wrong.chanlocs, 'maplimits', [-max(abs(ERP_wrong(:))), max(abs(ERP_wrong(:)))]);
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end

