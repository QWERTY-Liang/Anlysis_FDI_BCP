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

%load MuBeta_EOL_HL.mat %
%load Beta_EOL_HL.mat %
%load Mu_EOL_HL.mat %
%load Delta_EOL_HL.mat
%load Theta_EOL_HL.mat

%load MuBeta_EOL_FDIBCP.mat
%load Beta_EOL_FDIBCP.mat
%load Mu_EOL_FDIBCP.mat
%load Delta_EOL_FDIBCP.mat

%load MuBeta_RL_HL.mat
%load Beta_RL_HL.mat
%load Mu_RL_HL.mat
load Delta_RL_HL.mat

%load MuBeta_RL_FDIBCP.mat
%load Beta_RL_FDIBCP.mat
%load Mu_RL_FDIBCP.mat
%load Delta_RL_FDIBCP.mat

%fftlen = round(fs/21.5*6); % Window of how many sample points? If there is an SSVEP involved, whether or not you are interested in analyzing it, it is good to have all power related to the SSVEP isolated in a single frequency bin. This happens when you choose a window length that is an integer number of SSVEP cycles.
%F = [0:fftlen-1]*fs/fftlen; % frequency scale, given window length (remember resolution = 1/window-duration)
%ff = find((F>13 & F<21) | (F>22 & F<30)); % the indices of F that cover the spectral range/band of interest. Let's say we're interested in Mu and Beta bands combined. Note here I'm avoiding the SSVEP frequency (18.75hz in this example)
 Ts = [-2700:20:800]; % in msec, centered on what times do you want to measure spectral amplitude? i.e., where to center each consecutive window in time
% Tr = [-188*3:47:188]; % for response-locked
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
%全局看
num_plots = 15;
time_points = linspace(-188*3,188*1, num_plots);  % in milliseconds

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

%% Mu/beta waveforms 单频道检查
% %ch = [124 51];
% ch = [115 54];%[6 35]%[116 55]; % select left/right channels - typically D19 and B22, and that's the case here
% %ch = [124 50];%后 CP34
% %ch = [108 63];%前 FC34
% %ch = [116 55];%外 C34+
% %ch = [114 53];%里 C34-
% 
% figure; hold on
% 
% %correct+
% h1=plot(Ts,(mean(avMB1(ch(2),:,:),3)+mean(avMB2(ch(1),:,:),3))/2,'r'); % contralateral to correct side
% h2=plot(Ts,(mean(avMB1(ch(1),:,:),3)+mean(avMB2(ch(2),:,:),3))/2,'--r');%(mean(avMB(ch(1),:,1,c,sbj),5)+mean(avMB(ch(2),:,2,c,sbj),5))/2,'--','Color',colours(c,:),'LineWidth',2) % ipsilateral (dashed)
% %wrong
% h3=plot(Ts,(mean(avMB11(ch(2),:,:),3)+mean(avMB22(ch(1),:,:),3))/2,'b'); % contralateral to correct side
% h4=plot(Ts,(mean(avMB11(ch(1),:,:),3)+mean(avMB22(ch(2),:,:),3))/2,'--b');%
% set(gca,'Ydir','reverse') % we often turn the y axis upside down, to show increasing motor preparation (which is reflected in decreasing MB amplitude)
% %ylim([6.5,8.5])
% xlim([-500,200])
% title('MB: contralateral vs ipsilateral (dash)');
% %legend([h1,h2,h3,h4],{'Correct_contra','Correct_ipsi','Wrong_contra','Wrong_ipsi'},'AutoUpdate', 'off');
% %legend([h1,h2,h3,h4],{'High contra','High ipsi','Low contra','Low ipsi'},'AutoUpdate', 'off');
% legend([h1,h2,h3,h4],{'High contra','High ipsi','Low contra','Low ipsi'},'AutoUpdate', 'off');
% %legend([h1,h2,h3,h4],{'FDI contra','FDI ipsi','BCP contra','BCP ipsi'},'AutoUpdate', 'off');
% 
% xlabel('Time (ms)');
% ylabel('Beta (µV/m^2)');
% title(['MB for Channel :', num2str(ch(1)),' and : ', num2str(ch(2))]);
% 
% % Add vertical lines at specified time points
% xline([0, 100], '--r', { 'EMG onsite', 'appx. respT'});
%%
% 选择相邻电

% % 左侧（115）和右侧（54）附近电极
% left_adjacent = [114 115 116 108 124];  % 示例相邻电极
% right_adjacent = [53 54 55 63 50];     % 示例相邻电极
% 
% figure; hold on
% 
% % 计算并绘制 correct 反应
% % contralateral to correct side
% h1 = plot(Ts, (mean(mean(avMB1(right_adjacent,:,:),1),3) + mean(mean(avMB2(left_adjacent,:,:),1),3)) / 2, 'r'); 
% % ipsilateral (dashed)
% h2 = plot(Ts, (mean(mean(avMB1(left_adjacent,:,:),1),3) + mean(mean(avMB2(right_adjacent,:,:),1),3)) / 2, '--r');
% 
% % 计算并绘制 wrong 反应
% % contralateral to correct side
% h3 = plot(Ts, (mean(mean(avMB11(right_adjacent,:,:),1),3) + mean(mean(avMB22(left_adjacent,:,:),1),3)) / 2, 'b');
% % ipsilateral (dashed)
% h4 = plot(Ts, (mean(mean(avMB11(left_adjacent,:,:),1),3) + mean(mean(avMB22(right_adjacent,:,:),1),3)) / 2, '--b');
% set(gca,'Ydir','reverse') % we often turn the y axis upside down, to show increasing motor preparation (which is reflected in decreasing MB amplitude)
% %ylim([6.5,8.5])
% xlim([-600,200])
% title('MB: contralateral vs ipsilateral (dash)');
% %legend([h1,h2,h3,h4],{'Correct_contra','Correct_ipsi','Wrong_contra','Wrong_ipsi'},'AutoUpdate', 'off');
% %legend([h1,h2,h3,h4],{'High contra','High ipsi','Low contra','Low ipsi'},'AutoUpdate', 'off');
% %legend([h1,h2,h3,h4],{'High contra','High ipsi','Low contra','Low ipsi'},'AutoUpdate', 'off');
% legend([h1,h2,h3,h4],{'FDI contra','FDI ipsi','BCP contra','BCP ipsi'},'AutoUpdate', 'off');
% 
% xlabel('Time (ms)');
% ylabel('Beta (µV/m^2)');
% title(['MB for Channel avg']);
% % Add vertical lines at specified time points
% xline([0, 100], '--r', { 'EMG onsite', 'appx. respT'});
%%
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
h=TLBF3_stdshade((mean(avMB1(right_adjacent,:,:), 3) + mean(avMB2(left_adjacent,:,:), 3)) / 2, 0.15, color_correct1, Ts, smooth);

% Ipsilateral (dashed)
h2=TLBF3_stdshade_dash((mean(avMB1(left_adjacent,:,:), 3) + mean(avMB2(right_adjacent,:,:), 3)) / 2, 0.15, color_correct2, Ts, smooth);

% Calculate and plot wrong responses with shading
% Contralateral to correct side
h3=TLBF3_stdshade((mean(avMB11(right_adjacent,:,:), 3) + mean(avMB22(left_adjacent,:,:), 3)) / 2, 0.15, color_wrong1, Ts, smooth);

% Ipsilateral (dashed)
h4=TLBF3_stdshade_dash((mean(avMB11(left_adjacent,:,:), 3) + mean(avMB22(right_adjacent,:,:), 3)) / 2, 0.15, color_wrong2, Ts, smooth);

set(gca, 'Ydir', 'reverse');  % Reverse the Y-axis

title('MB: contralateral vs ipsilateral (dash)');
 %legend({'FDI contra', 'FDI ipsi', 'BCP contra', 'BCP ipsi'}, 'AutoUpdate', 'off');
legend({'High contra', 'High ipsi', 'Low contra', 'Low ipsi'}, 'AutoUpdate', 'off');

xlabel('Time (ms)');
ylabel('MuBeta (µV/m^2)');

% Add vertical lines at specified time points
%xline([0, 100], '--r', {'EMG onset', 'appx. respT'});xlim([-600, 200]);
xline([-120 0], '--r', {'appx.rsepT', 'EVend'});xlim([-600, 200]);%ylim([10,12]);

%% % Plot the scalp topographies at selected time points
%选平均时间段 -20ms-280ms
%Ts_idx=135:150;%FDIBCP RL
figure
%Ts_idx=125:141;%HL EOL
Ts_idx=120:137; % FDIBCP EOL

    subplot(1,2, 1);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB1(1:128, Ts_idx,:),3),2) - mean(mean(avMB2(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs,'electrodes', 'numbers');%,'electrodes', 'labels');
    title('FDI');
    %title('High');
    colorbar;

    subplot(1,2, 2);  % Create a 2x5 subplot
    topoplot( double(  mean(mean(avMB11(1:128, Ts_idx,:),3),2) - mean(mean(avMB22(1:128, Ts_idx,:),3),2)   ), ...
        EEG.chanlocs,'electrodes', 'numbers');%,'electrodes', 'labels');
     title('BCP');
    %title('Low');
    colorbar;

