%%
%1.每人单独分析
%2. 加了'调整allbehaviour'以及检查
%% Toolbox requirements: 
clc
clear
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2'));% current folder

% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);
eeglab
%
%%名称规则 波段_切分方法_分类方法_肌肉/条件
%filepath='G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\Motor MuBeta\RT_split3';
cd ('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\Motor MuBeta\RT_split3\');
%load MuBeta_RL_HL.mat  %
%load Delta_RL_HL.mat  %
%load MuBeta_RL_FDIBCP.mat  %

%load Delta_RL_FDIBCP.mat  %

%load MuBeta_EOL_HL.mat  %
%load Delta_EOL_HL.mat  %
%load MuBeta_EOL_FDIBCP.mat

load Delta_EOL_FDIBCP.mat

Ts = [-2700:20:800];%RL和EoL时的时间坐标

load 'chanlocsBioSemi_128_EOG4_Liang.mat'
EEG.chanlocs=chanlocs;

%% % Plot the scalp topographies at selected time points
%全局看 快速RT

avMB1=avMB1_1;
avMB2=avMB2_1;
avMB11=avMB11_1;
avMB22=avMB22_1;

num_plots = 15;
% RT的分段
time_points = linspace(-188*3,188*1, num_plots);  % in milliseconds
%SL的分段
%time_points = linspace(-900,1600, num_plots);
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
time_points = linspace(-188*3,188*1, num_plots);  % in milliseconds
%SL的分段
% time_points = linspace(-900,1600, num_plots);
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
time_points = linspace(-188*3,188*1, num_plots);  % in milliseconds
%SL的分段
% time_points = linspace(-900,1600, num_plots);
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
%%
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

%沙滩色
% #96ceb4	(150,206,180)
% #ffeead	(255,238,173)
% #ff6f69	(255,111,105)
% #ffcc5c	(255,204,92)
% #88d8b0	(136,216,176)
left_adjacent = [114 115 116 108 124 ];  % 示例相邻电极
right_adjacent = [53 54 55 63 50 ];     % 示例相邻电极

color_correct1 = [255,111,105]/255;  % Red for correct responses

color_correct1_1 = [255,204,92]/255;  % Red for correct responses

color_correct1_2 = [136,216,176]/255;  % Red for correct responses


   % Blue for wrong responses
color_wrong1 = [255,111,105]/255;
 % Blue for wrong responses
color_wrong1_1 = [255,204,92]/255;
    % Blue for wrong responses
color_wrong1_2 = [136,216,176]/255;

figure('units','normalized','outerposition',[0 0 1 1]); 
subplot(1,2,1)

hold on
smooth=4;%这里改变平滑曲线
% Calculate and plot correct responses with shading (stdshade function needed)
% Contralateral to correct side
h=TLBF3_stdshade( squeeze(  (mean(avMB1_1(right_adjacent,:,:), 1) + mean(avMB2_1(left_adjacent,:,:), 1)) / 2)' , 0.15, color_correct1, Ts, smooth);
h_1=TLBF3_stdshade(squeeze((mean(avMB1_2(right_adjacent,:,:), 1) + mean(avMB2_2(left_adjacent,:,:), 1)) / 2)', 0.15, color_correct1_1, Ts, smooth);
%h_2=TLBF3_stdshade_dash(squeeze((mean(avMB1_3(right_adjacent,:,:), 1) + mean(avMB2_3(left_adjacent,:,:), 1)) / 2)', 0.15, color_correct1_2, Ts, smooth);
h_2=TLBF3_stdshade(squeeze((mean(avMB1_3(right_adjacent,:,:), 1) + mean(avMB2_3(left_adjacent,:,:), 1)) / 2)', 0.15, color_correct1_2, Ts, smooth);
xlabel('Time (ms)', 'FontSize',20);
ylabel('MuBeta (µV/m^2)', 'FontSize',20);
legend({'Q-RT', 'M-RT', 'S-RT'}, 'AutoUpdate', 'off', 'FontSize',15);
set(gca, 'Ydir', 'reverse', 'FontSize',20);  % Reverse the Y-axis
% Add vertical lines at specified time points
%xline([-232, 0], '--r', {'appx.EMG onset', 'RL'});xlim([-1000, 300]);%ylim([10.5,12.5]);
xline([0 100], '--r', {'EMG onset', 'appx.rsepT'});xlim([-900, 400]);%ylim([10,12]);
%title('High contrast RT- split3 (contralateral)', 'FontSize',20);
title('FDI muscle RT- split3 (contralateral)', 'FontSize',20);
hold off


subplot(1,2,2)
hold on
%h3=TLBF3_stdshade_dash(squeeze((mean(avMB11_1(right_adjacent,:,:), 1) + mean(avMB22_1(left_adjacent,:,:), 1)) / 2)', 0.15, color_wrong1, Ts, smooth);
h3=TLBF3_stdshade(squeeze((mean(avMB11_1(right_adjacent,:,:), 1) + mean(avMB22_1(left_adjacent,:,:), 1)) / 2)', 0.15, color_wrong1, Ts, smooth);
h3_1=TLBF3_stdshade(squeeze((mean(avMB11_2(right_adjacent,:,:), 1) + mean(avMB22_2(left_adjacent,:,:), 1)) / 2)', 0.15, color_wrong1_1, Ts, smooth);
h3_2=TLBF3_stdshade(squeeze((mean(avMB11_3(right_adjacent,:,:), 1) + mean(avMB22_3(left_adjacent,:,:), 1)) / 2)', 0.15, color_wrong1_2, Ts, smooth);

xlabel('Time (ms)', 'FontSize',20);
ylabel('MuBeta (µV/m^2)', 'FontSize',20);
legend({'Q-RT', 'M-RT', 'S-RT'}, 'AutoUpdate', 'off', 'FontSize',15);
set(gca, 'Ydir', 'reverse', 'FontSize',20);  % Reverse the Y-axis
% Add vertical lines at specified time points
%xline([-232, 0], '--r', {'appx.EMG onset', 'RL'});xlim([-1000, 300]);%ylim([10.5,12.5]);
xline([0 100], '--r', {'EMG onset', 'appx.rsepT'});xlim([-900, 400]);%ylim([10,12]);
%title('Low contrast RT- split3 (contralateral)', 'FontSize',20);
title('BCP muscle RT- split3 (contralateral)', 'FontSize',20);
hold off


%% % Plot the scalp topographies at selected time points
%选平均时间段 -20ms-280ms
%Ts_idx=135:150;%FDIBCP RL

avMB1=avMB1_1;
avMB2=avMB2_1;
avMB11=avMB11_1;
avMB22=avMB22_1;

% min_val = -0.6;
% max_val =0.6;

t1='FDI';
t2='BCP';

figure
Ts_idx=116:137;%HL EOL


    subplot(3,2, 1);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB1(1:128, Ts_idx,:),3),2) - mean(mean(avMB2(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    %title('FDI (quick RT)');
     title([t1 '(Q-RT)'], 'FontSize',12);
     % Set the colormap (e.g., 'jet', 'hot', etc.)
 %colormap(cool); % You can change 'jet' to your preferred colormap
% 
% % Set the color limits (adjust min_val and max_val based on your data range)
% caxis([min_val max_val]);

% Display the colorbar
colorbar;


    subplot(3,2, 2);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB11(1:128, Ts_idx,:),3),2) - mean(mean(avMB22(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    % title('BCP (quick RT)');
     title([t2 '(Q-RT)'], 'FontSize',12);
    colorbar;
% Plot the scalp topographies at selected time points
%选平均时间段 -20ms-280ms
%Ts_idx=135:150;%FDIBCP RL

avMB1=avMB1_2;
avMB2=avMB2_2;
avMB11=avMB11_2;
avMB22=avMB22_2;
Ts_idx=116:137;%HL EOL


    subplot(3,2, 3);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB1(1:128, Ts_idx,:),3),2) - mean(mean(avMB2(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    %title('FDI');
     title([t1 '(M-RT)'], 'FontSize',12);
    colorbar;

    subplot(3,2, 4);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB11(1:128, Ts_idx,:),3),2) - mean(mean(avMB22(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    % title('BCP');
     title([t1 '(M-RT)'], 'FontSize',12);
    colorbar;
Ts_idx=116:137;%HL EOL

avMB1=avMB1_3;
avMB2=avMB2_3;
avMB11=avMB11_3;
avMB22=avMB22_3;

%Ts_idx=91:128; % quick RT 880-1600ms for Delta SL
%Ts_idx=89:120; % quick RT 800-1400ms/ 350-530 for Delta SL FDIBCP
%Ts_idx=66:75; % quick RT 800-1400ms/ 350-530 for Delta SL FDIBCP %BCP时间提前且固定
%figure
%Ts_idx=125:141;%HL EOL


    subplot(3,2, 5);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB1(1:128, Ts_idx,:),3),2) - mean(mean(avMB2(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    %title('FDI');
     title([t1 '(S-RT)'], 'FontSize',12);
    colorbar;

    subplot(3,2, 6);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB11(1:128, Ts_idx,:),3),2) - mean(mean(avMB22(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    % title('BCP');
     title([t1 '(S-RT)'], 'FontSize',12);
    colorbar;


    %% 补充总图 3个RT一起
%测试stdshade
% Define colors for shading
%蓝色浅到深
% #011f4b	(1,31,75)
% #03396c	(3,57,108)
% #005b96	(0,91,150)
% #6497b1	(100,151,177)
% #b3cde0	(179,205,224)
%粉色浅到深
% #ffc2cd	(255,194,205)
% #ff93ac	(255,147,172)
% #ff6289	(255,98,137)
% #fc3468	(252,52,104)
% #ff084a	(255,8,74)
left_adjacent = [114 115 116 108 124];  % 示例相邻电极
right_adjacent = [53 54 55 63 50];     % 示例相邻电极
color_correct1 = [3,57,108]/255;  % Red for correct responses
color_correct2 = [100,151,177]/255; 
color_wrong1 = [255,8,74]/255;    % Blue for wrong responses
color_wrong2 = [255,194,205]/255;
figure; hold on
smooth=4;%这里改变平滑曲线
% Calculate and plot correct responses with shading (stdshade function needed)
% Contralateral to correct side

h=TLBF3_stdshade(squeeze((    mean(avMB1_1(right_adjacent,:,:), 1) + mean(avMB2_1(left_adjacent,:,:), 1) +  mean(avMB1_2(right_adjacent,:,:), 1) + mean(avMB2_2(left_adjacent,:,:), 1) +mean(avMB1_3(right_adjacent,:,:), 1) + mean(avMB2_3(left_adjacent,:,:), 1)  ) / (2*3))', 0.15, color_correct1, Ts, smooth);

% Ipsilateral (dashed)
h2=TLBF3_stdshade_dash( squeeze(   (   mean(avMB1_1(left_adjacent,:,:), 1) + mean(avMB2_1(right_adjacent,:,:), 1) +  mean(avMB1_2(left_adjacent,:,:), 1) + mean(avMB2_2(right_adjacent,:,:), 1) +  mean(avMB1_3(left_adjacent,:,:), 1) + mean(avMB2_3(right_adjacent,:,:), 1)     ) / (2*3))', 0.15, color_correct2, Ts, smooth);

% Calculate and plot wrong responses with shading
% Contralateral to correct side
h3=TLBF3_stdshade(   squeeze(   (   mean(avMB11_1(right_adjacent,:,:), 1) + mean(avMB22_1(left_adjacent,:,:), 1)  +  mean(avMB11_2(right_adjacent,:,:), 1) + mean(avMB22_2(left_adjacent,:,:), 1)   +  mean(avMB11_3(right_adjacent,:,:), 1) + mean(avMB22_3(left_adjacent,:,:), 1)    ) / (2*3))', 0.15, color_wrong1, Ts, smooth);

% Ipsilateral (dashed)
h4=TLBF3_stdshade_dash(  squeeze(      (   mean(avMB11_1(left_adjacent,:,:), 1) + mean(avMB22_1(right_adjacent,:,:), 1) +    mean(avMB11_2(left_adjacent,:,:), 1) + mean(avMB22_2(right_adjacent,:,:), 1)   +    mean(avMB11_3(left_adjacent,:,:), 1) + mean(avMB22_3(right_adjacent,:,:), 1) ) / (2*3))'  , 0.15, color_wrong2, Ts, smooth);

set(gca, 'Ydir', 'reverse');  % Reverse the Y-axis

title('MB: contralateral vs ipsilateral (dash)');
 legend({'FDI contra', 'FDI ipsi', 'BCP contra', 'BCP ipsi'}, 'AutoUpdate', 'off');
%legend({'High contra', 'High ipsi', 'Low contra', 'Low ipsi'}, 'AutoUpdate', 'off');

xlabel('Time (ms)');
ylabel('MuBeta (µV/m^2)');

% Add vertical lines at specified time points
xline([0, 100], '--r', {'EMG onset', 'appx. respT'});xlim([-600, 200]);
%xline([-120 0], '--r', {'appx.rsepT', 'EVend'});xlim([-600, 200]);%ylim([10,12]);