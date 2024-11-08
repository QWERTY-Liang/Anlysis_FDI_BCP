%%
%1.每人单独分析
%2. 加了'调整allbehaviour'以及检查
%% Toolbox requirements: 
clc
clear
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2'));% current folder

%% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);
eeglab
%%
%%名称规则 波段_切分方法_分类方法_肌肉/条件
%filepath='G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\Motor MuBeta\RT_split3';
cd ('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\Motor MuBeta\RT_split3\');
 load MuBeta_SL_HL.mat  %
%load Delta_SL_HL.mat  %
%load MuBeta_SL_FDIBCP.mat  %
%load Delta_SL_FDIBCP.mat  %

Ts = [-950:20:2700]; % SL时的时间坐标


%Ts = [-2700:20:800];RL和EoL时的时间坐标

%fftlen = round(fs/21.5*6); % Window of how many sample points? If there is an SSVEP involved, whether or not you are interested in analyzing it, it is good to have all power related to the SSVEP isolated in a single frequency bin. This happens when you choose a window length that is an integer number of SSVEP cycles.
%F = [0:fftlen-1]*fs/fftlen; % frequency scale, given window length (remember resolution = 1/window-duration)
%ff = find((F>13 & F<21) | (F>22 & F<30)); % the indices of F that cover the spectral range/band of interest. Let's say we're interested in Mu and Beta bands combined. Note here I'm avoiding the SSVEP frequency (18.75hz in this example)
  % in msec, centered on what times do you want to measure spectral amplitude? i.e., where to center each consecutive window in time


load 'chanlocsBioSemi_128_EOG4_Liang.mat'
EEG.chanlocs=chanlocs;
       %% Now let's look at Mu/Beta 'MB'

% Let's first verify that there is lateralisation of MB, where the amplitude is lower contralateral to the button that was pressed. Plotting
% this topography also serves to highlight which electrodes might be best for measuring MB (although note there are some tasks like the
% continuous dots (2013 paper) and any of our delayed-response tasks, where there is precious little difference in MB amplitude contra/ipsi
% to responding hand, i.e. not much lateralisation):
% Define the number of plots and the time points
% sub=1;%每个人检查
% num_plots = 15;
% time_points = linspace(-1600,900, num_plots);  % in milliseconds
% 
% % Find the indices corresponding to the defined time points
% time_indices = arrayfun(@(t) find(Ts >= t, 1), time_points);
% %只画左和右, beta层面baseline 8成功
% figure;
% for i = 1:num_plots
%     subplot(3, 5, i);  % Create a 2x5 subplot
%     topoplot(   double(avMB1(1:128, time_indices(i),sub)-avMB2(1:128, time_indices(i),sub))    , ...
%         EEG.chanlocs,  'colormap','jet');%,'electrodes', 'labels');%
%     title([num2str(time_points(i)), ' ms']);
%     colorbar;
% end
% 
% figure;
% for i = 1:num_plots
%     subplot(3, 5, i);  % Create a 2x5 subplot
%     topoplot(   double(avMB11(1:128, time_indices(i),sub)-avMB22(1:128, time_indices(i),sub))    , ...
%         EEG.chanlocs,  'colormap','jet');%,'electrodes', 'labels');%
%     title([num2str(time_points(i)), ' ms']);
%     colorbar;
% end


%% % Plot the scalp topographies at selected time points
%全局看 快速RT
avMB1=avMB1_1;
avMB2=avMB2_1;
avMB11=avMB11_1;
avMB22=avMB22_1;

num_plots = 15;
% RT的分段
%time_points = linspace(-188*3,188*1, num_plots);  % in milliseconds
%SL的分段
time_points = linspace(-900,1600, num_plots);
% Find the indices corresponding to the defined time points
time_indices = arrayfun(@(t) find(Ts >= t, 1), time_points);
figure;
for i = 1:num_plots
    subplot(3, 5, i);  % Create a 2x5 subplot
    topoplot( double(  mean(avMB1(1:128, time_indices(i),:),3) - mean(avMB2(1:128, time_indices(i),:),3)   ), ...
        EEG.chanlocs,'electrodes', 'numbers');%,'electrodes', 'labels');
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end


time_indices = arrayfun(@(t) find(Ts >= t, 1), time_points);
figure;
for i = 1:num_plots
    subplot(3, 5, i);  % Create a 2x5 subplot
    topoplot( double(  mean(avMB11(1:128, time_indices(i),:),3) - mean(avMB22(1:128, time_indices(i),:),3)   ), ...
        EEG.chanlocs,'electrodes', 'numbers');%,'electrodes', 'labels');
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end
% trange = find(Ts>=800,1); % pick a time just before response
% figure;
% topoplot(double(avMB1(1:128,trange)-avMB2(1:128,trange)),EEG.chanlocs,'electrodes','labels','colormap','jet');%,'maplimits',scale));
% %topoplot(double(mean(mean(mean(avMB(1:nchan,trange,1,:,sbj),5),4),2)-mean(mean(mean(avMBr(1:nchan,trange,2,:,sbj),5),4),2)),chanlocs,'electrodes','labels','colormap','jet');%,'maplimits',scale) % some additional arguments that you can add: ,'maplimits',[-1 1]*4 sets the colorscale
% %topoplot(ERP_right(:, time_indices(i)), EEG_selected_right.chanlocs, 'maplimits', [-max(abs(ERP_right(:))), max(abs(ERP_right(:)))]);
% 
% title('Left minus Right press')
% colorbar
% Plot the scalp topographies at selected time points

%% % Plot the scalp topographies at selected time points
%全局看 中速RT
avMB1=avMB1_2;
avMB2=avMB2_2;
avMB11=avMB11_2;
avMB22=avMB22_2;

num_plots = 15;
% RT的分段
%time_points = linspace(-188*3,188*1, num_plots);  % in milliseconds
%SL的分段
time_points = linspace(-900,1600, num_plots);
% Find the indices corresponding to the defined time points
time_indices = arrayfun(@(t) find(Ts >= t, 1), time_points);
figure;
for i = 1:num_plots
    subplot(3, 5, i);  % Create a 2x5 subplot
    topoplot( double(  mean(avMB1(1:128, time_indices(i),:),3) - mean(avMB2(1:128, time_indices(i),:),3)   ), ...
        EEG.chanlocs,'electrodes', 'numbers');%,'electrodes', 'labels');
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end


time_indices = arrayfun(@(t) find(Ts >= t, 1), time_points);
figure;
for i = 1:num_plots
    subplot(3, 5, i);  % Create a 2x5 subplot
    topoplot( double(  mean(avMB11(1:128, time_indices(i),:),3) - mean(avMB22(1:128, time_indices(i),:),3)   ), ...
        EEG.chanlocs,'electrodes', 'numbers');%,'electrodes', 'labels');
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end
%% % Plot the scalp topographies at selected time points
%全局看 慢速RT
avMB1=avMB1_3;
avMB2=avMB2_3;
avMB11=avMB11_3;
avMB22=avMB22_3;

num_plots = 15;
% RT的分段
%time_points = linspace(-188*3,188*1, num_plots);  % in milliseconds
%SL的分段
time_points = linspace(-900,1600, num_plots);
% Find the indices corresponding to the defined time points
time_indices = arrayfun(@(t) find(Ts >= t, 1), time_points);
figure;
for i = 1:num_plots
    subplot(3, 5, i);  % Create a 2x5 subplot
    topoplot( double(  mean(avMB1(1:128, time_indices(i),:),3) - mean(avMB2(1:128, time_indices(i),:),3)   ), ...
        EEG.chanlocs,'electrodes', 'numbers');%,'electrodes', 'labels');
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end


time_indices = arrayfun(@(t) find(Ts >= t, 1), time_points);
figure;
for i = 1:num_plots
    subplot(3, 5, i);  % Create a 2x5 subplot
    topoplot( double(  mean(avMB11(1:128, time_indices(i),:),3) - mean(avMB22(1:128, time_indices(i),:),3)   ), ...
        EEG.chanlocs,'electrodes', 'numbers');%,'electrodes', 'labels');
    title([num2str(time_points(i)), ' ms']);
    colorbar;
end
%% RT split3
%测试stdshade
% Define colors for shading
%蓝色浅到深
% #011f4b	(1,31,75)
% #03396c	(3,57,108)
% #005b96	(0,91,150)
% #6497b1	(100,151,177)
% #b3cde0	(179,205,224)
%粉色深到浅
% #ffc2cd	(255,194,205)
% #ff93ac	(255,147,172)
% #ff6289	(255,98,137)
% #fc3468	(252,52,104)
% #ff084a	(255,8,74)
left_adjacent = [114 115 116 108 124];  % 示例相邻电极
right_adjacent = [53 54 55 63 50];     % 示例相邻电极

color_correct1 = [1,31,75]/255;  % Red for correct responses
color_correct2 = [0,91,150]/255; 
color_correct1_1 = [3,57,108]/255;  % Red for correct responses
color_correct2_1 = [100,151,177]/255; 
color_correct1_2 = [0,91,150]/255;  % Red for correct responses
color_correct2_2 = [179,205,224]/255;

color_wrong2 = [255,194,205]/255;    % Blue for wrong responses
color_wrong1 = [255,98,137]/255;
color_wrong2_1 = [255,147,172]/255;    % Blue for wrong responses
color_wrong1_1 = [252,52,104]/255;
color_wrong2_2 = [255,98,137]/255;    % Blue for wrong responses
color_wrong1_2 = [255,8,74]/255;

figure(10); 
subplot(3,1,1)

hold on
smooth=4;%这里改变平滑曲线
% Calculate and plot correct responses with shading (stdshade function needed)
% Contralateral to correct side
h=TLBF3_stdshade((mean(avMB1_1(right_adjacent,:,:), 3) + mean(avMB2_1(left_adjacent,:,:), 3)) / 2, 0.15, color_correct1, Ts, smooth);

% Ipsilateral (dashed)
h2=TLBF3_stdshade_dash((mean(avMB1_1(left_adjacent,:,:), 3) + mean(avMB2_1(right_adjacent,:,:), 3)) / 2, 0.15, color_correct2, Ts, smooth);
% Calculate and plot wrong responses with shading
% Contralateral to correct side
h3=TLBF3_stdshade((mean(avMB11_1(right_adjacent,:,:), 3) + mean(avMB22_1(left_adjacent,:,:), 3)) / 2, 0.15, color_wrong1, Ts, smooth);

% Ipsilateral (dashed)
h4=TLBF3_stdshade_dash((mean(avMB11_1(left_adjacent,:,:), 3) + mean(avMB22_1(right_adjacent,:,:), 3)) / 2, 0.15, color_wrong2, Ts, smooth);

set(gca, 'Ydir', 'reverse');  % Reverse the Y-axis
%xlim([-800, 400]);
title('Quick RT: contralateral vs ipsilateral (dash)');
 legend({'FDI contra', 'FDI ipsi', 'BCP contra', 'BCP ipsi'}, 'AutoUpdate', 'off');
% legend({'High contra', 'High ipsi', 'Low contra', 'Low ipsi'}, 'AutoUpdate', 'off');

xlabel('Time (ms)');
ylabel('MuBeta (µV/m^2)');
% Add vertical lines at specified time points
xline([-600, 0, 800], '--r', {'cue on','evidence on', 'minEvd0.8'});xlim([-650, 1500]);
hold off
%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2)

hold on

h_1=TLBF3_stdshade((mean(avMB1_2(right_adjacent,:,:), 3) + mean(avMB2_2(left_adjacent,:,:), 3)) / 2, 0.15, color_correct1_1, Ts, smooth);

% Ipsilateral (dashed)
h2_1=TLBF3_stdshade_dash((mean(avMB1_2(left_adjacent,:,:), 3) + mean(avMB2_2(right_adjacent,:,:), 3)) / 2, 0.15, color_correct2_1, Ts, smooth);
h3_1=TLBF3_stdshade((mean(avMB11_2(right_adjacent,:,:), 3) + mean(avMB22_2(left_adjacent,:,:), 3)) / 2, 0.15, color_wrong1_1, Ts, smooth);

% Ipsilateral (dashed)
h4_1=TLBF3_stdshade_dash((mean(avMB11_2(left_adjacent,:,:), 3) + mean(avMB22_2(right_adjacent,:,:), 3)) / 2, 0.15, color_wrong2_1, Ts, smooth);

set(gca, 'Ydir', 'reverse');  % Reverse the Y-axis
%xlim([-800, 400]);
title('Mid RT: contralateral vs ipsilateral (dash)');
 legend({'FDI contra', 'FDI ipsi', 'BCP contra', 'BCP ipsi'}, 'AutoUpdate', 'off');
% legend({'High contra', 'High ipsi', 'Low contra', 'Low ipsi'}, 'AutoUpdate', 'off');

xlabel('Time (ms)');
ylabel('MuBeta (µV/m^2)');
% Add vertical lines at specified time points
xline([-600, 0, 800], '--r', {'cue on','evidence on', 'minEvd0.8'});xlim([-650, 1500]);
hold off


%%%%%%%%%%%%%%%%%
subplot(3,1,3)

hold on

h_2=TLBF3_stdshade((mean(avMB1_3(right_adjacent,:,:), 3) + mean(avMB2_3(left_adjacent,:,:), 3)) / 2, 0.15, color_correct1_2, Ts, smooth);

% Ipsilateral (dashed)
h2_2=TLBF3_stdshade_dash((mean(avMB1_3(left_adjacent,:,:), 3) + mean(avMB2_3(right_adjacent,:,:), 3)) / 2, 0.15, color_correct2_2, Ts, smooth);
h3_2=TLBF3_stdshade((mean(avMB11_3(right_adjacent,:,:), 3) + mean(avMB22_3(left_adjacent,:,:), 3)) / 2, 0.15, color_wrong1_2, Ts, smooth);

% Ipsilateral (dashed)
h4_2=TLBF3_stdshade_dash((mean(avMB11_3(left_adjacent,:,:), 3) + mean(avMB22_3(right_adjacent,:,:), 3)) / 2, 0.15, color_wrong2_2, Ts, smooth);


set(gca, 'Ydir', 'reverse');  % Reverse the Y-axis
%xlim([-800, 400]);
title('Slow RT: contralateral vs ipsilateral (dash)');
 legend({'FDI contra', 'FDI ipsi', 'BCP contra', 'BCP ipsi'}, 'AutoUpdate', 'off');
% legend({'High contra', 'High ipsi', 'Low contra', 'Low ipsi'}, 'AutoUpdate', 'off');

xlabel('Time (ms)');
ylabel('MuBeta (µV/m^2)');
% Add vertical lines at specified time points
xline([-600, 0, 800], '--r', {'cue on','evidence on', 'minEvd0.8'});xlim([-650, 1500]);
hold off
% Add vertical lines at specified time points
%xline([0, 100], '--r', {'EMG onset', 'appx. respT'});

%% % Plot the scalp topographies at selected time points
%选平均时间段 -20ms-280ms
%Ts_idx=135:150;%FDIBCP RL

avMB1=avMB1_1;
avMB2=avMB2_1;
avMB11=avMB11_1;
avMB22=avMB22_1;
%Ts_idx=49:70; % quick RT 0-430ms for MB SL HL
%Ts_idx=57:91; % quick RT 170-850ms for Delta SL HL
Ts_idx=64:84; % quick RT 300-700ms for Delta SL FDIBCP



figure
%Ts_idx=125:141;%HL EOL


    subplot(1,2, 1);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB1(1:128, Ts_idx,:),3),2) - mean(mean(avMB2(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    title('FDI (quick RT)');
    % title('High (quick _ RT)');
    colorbar;

    subplot(1,2, 2);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB11(1:128, Ts_idx,:),3),2) - mean(mean(avMB22(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
     title('BCP (quick RT)');
    % title('Low (quick _ RT)');
    colorbar;
%% % Plot the scalp topographies at selected time points
%选平均时间段 -20ms-280ms
%Ts_idx=135:150;%FDIBCP RL

avMB1=avMB1_2;
avMB2=avMB2_2;
avMB11=avMB11_2;
avMB22=avMB22_2;
% Ts_idx=71:90; % Mid RT 600-900ms  SL HL
%Ts_idx=74:108; % quick RT 500-1200ms for Delta SL HL
Ts_idx=66:98; % quick RT 350-1000ms for Delta SL FDIBCP
figure
%Ts_idx=125:141;%HL EOL


    subplot(1,2, 1);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB1(1:128, Ts_idx,:),3),2) - mean(mean(avMB2(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    title('FDI');
    % title('High Mid _ RT');
    colorbar;

    subplot(1,2, 2);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB11(1:128, Ts_idx,:),3),2) - mean(mean(avMB22(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
     title('BCP');
    % title('Low Mid _ RT');
    colorbar;
%% % Plot the scalp topographies at selected time points
%选平均时间段 -20ms-280ms
%Ts_idx=135:150;%FDIBCP RL

avMB1=avMB1_3;
avMB2=avMB2_3;
avMB11=avMB11_3;
avMB22=avMB22_3;
% Ts_idx=90:102; % Slow RT 900-1200ms
%Ts_idx=91:128; % quick RT 880-1600ms for Delta SL
Ts_idx=89:120; % quick RT 800-1400ms/ 350-530 for Delta SL FDIBCP
%Ts_idx=66:75; % quick RT 800-1400ms/ 350-530 for Delta SL FDIBCP %BCP时间提前且固定
figure
%Ts_idx=125:141;%HL EOL


    subplot(1,2, 1);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB1(1:128, Ts_idx,:),3),2) - mean(mean(avMB2(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    title('FDI');
    % title('High Slow _ RT');
    colorbar;

    subplot(1,2, 2);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB11(1:128, Ts_idx,:),3),2) - mean(mean(avMB22(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
     title('BCP');
    % title('Low Slow _ RT');
    colorbar;