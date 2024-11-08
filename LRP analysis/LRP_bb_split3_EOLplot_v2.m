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

%cd ('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\Motor MuBeta\RT_split3\');


%load LRP_EOL_HL.mat
load LRP_EOL_FDIBCP.mat

load 'chanlocsBioSemi_128_EOG4_Liang.mat'
EEG.chanlocs=chanlocs;
%% plot for all RT
%沙滩色
% #96ceb4	(150,206,180)
% #ffeead	(255,238,173)
% #ff6f69	(255,111,105)
% #ffcc5c	(255,204,92)
% #88d8b0	(136,216,176)
% left_adjacent = [114 115 116 108 124 ];  % 示例相邻电极
% right_adjacent = [53 54 55 63 50 ];     % 示例相邻电极

color_correct1 = [255,111,105]/255;  % Red for correct responses

color_correct1_1 = [255,204,92]/255;  % Red for correct responses

color_correct1_2 = [136,216,176]/255;  % Red for correct responses


   % Blue for wrong responses
color_wrong1 = [255,111,105]/255;
 % Blue for wrong responses
color_wrong1_1 = [255,204,92]/255;
    % Blue for wrong responses
color_wrong1_2 = [136,216,176]/255;


figure;
smooth=4;%这里改变平滑曲线
% plot(time_vector, mean(erp1_1(channel_to_plot, :,:,:),4)   );
h=TLBF3_stdshade( squeeze(   mean(LRP_con1(:, :,:,:),4)   )' , 0.15, color_correct1, time_vector, smooth);
hold on;
% plot(time_vector, mean(mean(erp2_1(channel_to_plot, :,:,:),4),3)  );
h1=TLBF3_stdshade( squeeze(   mean(LRP_con2(:, :,:,:),4)   )' , 0.15, color_correct1_1, time_vector, smooth);


%legend('High(0.14)','Low(0.07)','AutoUpdate', 'off', 'FontSize',15)
legend('FDI','BCP','AutoUpdate', 'off', 'FontSize',15)
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(['LRP contra-lpsi (EOL)']);
hold off;
% Add vertical lines at specified time points
hold on;
xline([-1500, 0, 120], '--r', {'1.5 before EVoff','EMG onsite', 'appx. respT'});xlim([-600+120,220]);
hold off;

%% RT-split3 作画

figure('units','normalized','outerposition',[0 0 1 1]); 
subplot(1,2,1)
Ts=time_vector;
hold on
smooth=10;%这里改变平滑曲线
% Calculate and plot correct responses with shading (stdshade function needed)
% Contralateral to correct side
h=TLBF3_stdshade( squeeze(   LRP_con1(:, :,:,1) )' , 0.15, color_correct1, Ts, smooth);
h_1=TLBF3_stdshade(squeeze( LRP_con1(:, :,:,2) )', 0.15, color_correct1_1, Ts, smooth);
%h_2=TLBF3_stdshade_dash(squeeze((mean(avMB1_3(right_adjacent,:,:), 1) + mean(avMB2_3(left_adjacent,:,:), 1)) / 2)', 0.15, color_correct1_2, Ts, smooth);
h_2=TLBF3_stdshade(squeeze( LRP_con1(:, :,:,3) )', 0.15, color_correct1_2, Ts, smooth);
xlabel('Time (ms)', 'FontSize',20);
ylabel('LRP (µV)', 'FontSize',20);
legend({'Q-RT', 'M-RT', 'S-RT'}, 'AutoUpdate', 'off', 'FontSize',15);
%set(gca, 'Ydir', 'reverse', 'FontSize',20);  % Reverse the Y-axis
% Add vertical lines at specified time points
xline([-1500, 0, 120], '--r', {'1.5 before EVoff','EMG onsite', 'appx. respT'});xlim([-600+120,220]);
%title('High contrast RT- split3 (contralateral)', 'FontSize',20);
title('FDI muscle RT- split3 (contralateral)', 'FontSize',20);
hold off


subplot(1,2,2)
hold on
%h3=TLBF3_stdshade_dash(squeeze((mean(avMB11_1(right_adjacent,:,:), 1) + mean(avMB22_1(left_adjacent,:,:), 1)) / 2)', 0.15, color_wrong1, Ts, smooth);
h3=TLBF3_stdshade(squeeze(  LRP_con2(:, :,:,1)  )', 0.15, color_wrong1, Ts, smooth);
h3_1=TLBF3_stdshade(squeeze(  LRP_con2(:, :,:,2)   )', 0.15, color_wrong1_1, Ts, smooth);
h3_2=TLBF3_stdshade(squeeze(   LRP_con2(:, :,:,3)   )', 0.15, color_wrong1_2, Ts, smooth);

xlabel('Time (ms)', 'FontSize',20);
ylabel('LRP (µV)', 'FontSize',20);
legend({'Q-RT', 'M-RT', 'S-RT'}, 'AutoUpdate', 'off', 'FontSize',15);
%set(gca, 'Ydir', 'reverse', 'FontSize',20);  % Reverse the Y-axis
% Add vertical lines at specified time points
xline([-1500, 0, 120], '--r', {'1.5 before EVoff','EMG onsite', 'appx. respT'});xlim([-600+120,220]);
%title('Low contrast RT- split3 (contralateral)', 'FontSize',20);
title('BCP muscle RT- split3 (contralateral)', 'FontSize',20);
hold off

%% 画平均cpp

Ts_idx=1421:1793; % quick RT 0-500ms 


% t1='High';
% t2='Low';
t1='FDI';
t2='BCP';

figure
%Ts_idx=125:141;%HL EOL


    subplot(2,1, 1);  % Create a 2x5 subplot
    topoplot( double( squeeze(   mean(mean(mean(ERP_right(:, Ts_idx,:,:),4),3),2)-mean(mean(mean(ERP_left(:, Ts_idx,:,:),4),3),2)  )   ), ...
        EEG.chanlocs);%, 'maplimits', [-40,40]);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    %title('FDI (all RT R-L)');
     title([t1 '(all-RT R-L)'], 'FontSize',12);
    colorbar;

    subplot(2,1, 2);  % Create a 2x5 subplot
    topoplot( double( squeeze(   mean(mean(mean(ERP_right2(:, Ts_idx,:,:),4),3),2)-mean(mean(mean(ERP_left2(:, Ts_idx,:,:),4),3),2) )  ), ...
        EEG.chanlocs);%, 'maplimits', [-40,40]);%,'electrodes', 'numbers');%,'electrodes', 'labels');
    % title('BCP (all RT R-L)');
     title([t2 ' (all-RT R-L)'], 'FontSize',12);
    colorbar;