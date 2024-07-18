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

%% 1. Load Sl_BB pre-cue-188
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
%% 3.1 left vs right using C3 and C4
%For left-hand responses, C4 is the contralateral (right hemisphere) signal and C3 is the ipsilateral (left hemisphere) signal.
%For right-hand responses, C3 is the contralateral (left hemisphere) signal and C4 is the ipsilateral (right hemisphere) signal.
%LRP=Contra−Ipsi
% Define the selected trials vector
% channel_to_plot = 19;
left_channel = 108;%115;  % C34通道号 (115,54) %FC34 (108 63)
right_channel = 63;%54; 

selected_trials_1 = zeros(length(AllBehaviour_SL_bb),1);  % left C4-C3
selected_trials_2 = zeros(length(AllBehaviour_SL_bb),1);  % right C3-C4
count=0;
% switch%%%%%%%%%%calculate correct first
%     case 'CW_SL_b' %correct vs wrong; pre cue SL
        for i=1:length(AllBehaviour_SL_bb)
            %select trials
            if AllBehaviour_SL_bb(i,8)==1 && AllBehaviour_SL_bb(i,4)==1 %&& AllBehaviour_SL_bb(i,9)==1%&& AllBehaviour_SL_bb(i,3)==0.14
                selected_trials_1(i)=1;
            elseif AllBehaviour_SL_bb(i,8)==2 && AllBehaviour_SL_bb(i,4)==1 %&& AllBehaviour_SL_bb(i,9)==2%&& AllBehaviour_SL_bb(i,3)==0.14
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
EEG_selected_left = pop_select(EEG, 'trial', selected_idx_1);
% Compute the ERP
ERP_left = mean(EEG_selected_left.data, 3);  % averaging across the third dimension (trials)

LRP_left = (ERP_left(right_channel, :) - ERP_left(left_channel, :)) / 2;%这里先平均再相减

%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_left = EEG_selected_left.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% Extract selected trials---for wrong
selected_idx_right = find(selected_trials_2 == 1);
EEG_selected_right = pop_select(EEG, 'trial', selected_idx_right);
% Compute the ERP
ERP_right = mean(EEG_selected_right.data, 3);  % averaging across the third dimension (trials)
LRP_right = (ERP_right(left_channel, :) - ERP_right(right_channel, :)) / 2;%这里先平均再相减
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_right = EEG_selected_right.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

LRP_left_correct=LRP_left;
LRP_right_correct=LRP_right;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot scalp plot ----- left-right _correct
num_plots = 10;
time_points = linspace(0, 1692, num_plots);  % in milliseconds

% Find the indices of the defined time points
time_indices = arrayfun(@(t) find(EEG_selected_left.times >= t, 1), time_points);

% Plot scalp topographies
figure;
for i = 1:num_plots
    subplot(2, 5, i);  % create a 2x4 subplot
    topoplot(ERP_left(:, time_indices(i))-ERP_right(:, time_indices(i)), EEG_selected_left.chanlocs, 'maplimits', [-30,30]);
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end
%% find all wrong
       for i=1:length(AllBehaviour_SL_bb)
            %select trials
            if AllBehaviour_SL_bb(i,8)==1 && AllBehaviour_SL_bb(i,4)==2%&& AllBehaviour_SL_bb(i,9)==1%&& AllBehaviour_SL_bb(i,3)==0.14
                selected_trials_1(i)=1;
            elseif AllBehaviour_SL_bb(i,8)==2 && AllBehaviour_SL_bb(i,4)==2%&& AllBehaviour_SL_bb(i,9)==2%&& AllBehaviour_SL_bb(i,3)==0.14
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
EEG_selected_left = pop_select(EEG, 'trial', selected_idx_1);
% Compute the ERP
ERP_left = mean(EEG_selected_left.data, 3);  % averaging across the third dimension (trials)

LRP_left = (ERP_left(right_channel, :) - ERP_left(left_channel, :)) / 2;%这里先平均再相减

%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_left = EEG_selected_left.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% Extract selected trials---for wrong
selected_idx_right = find(selected_trials_2 == 1);
EEG_selected_right = pop_select(EEG, 'trial', selected_idx_right);
% Compute the ERP
ERP_right = mean(EEG_selected_right.data, 3);  % averaging across the third dimension (trials)
LRP_right = (ERP_right(left_channel, :) - ERP_right(right_channel, :)) / 2;%这里先平均再相减
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
time_vector_right = EEG_selected_right.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

%% plot scalp plot ----- left-right _wrong
num_plots = 10;
time_points = linspace(0, 1692, num_plots);  % in milliseconds

% Find the indices of the defined time points
time_indices = arrayfun(@(t) find(EEG_selected_left.times >= t, 1), time_points);

% Plot scalp topographies
figure;
for i = 1:num_plots
    subplot(2, 5, i);  % create a 2x4 subplot
    topoplot(ERP_left(:, time_indices(i))-ERP_right(:, time_indices(i)), EEG_selected_left.chanlocs, 'maplimits', [-30,30]);
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end

%% plot LRPs
figure;
h1=plot(time_vector_left, -(LRP_left_correct+LRP_right_correct)/2,'r');
hold on;
%h2=plot(time_vector_right, LRP_right_correct,'b');
h3=plot(time_vector_left, -(LRP_left+LRP_right)/2,'b');
%h4=plot(time_vector_right, LRP_right,'--k');

legend([h1,h3],{'FDI','BCP'},'AutoUpdate', 'off');
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(['LRP for Channel C3:', num2str(left_channel),' and C4: ', num2str(right_channel)]);
% Add vertical lines at specified time points
xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s'});
hold off;