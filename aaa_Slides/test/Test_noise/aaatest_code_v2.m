clc
clear all
%%
trial_num=6;

%test becips
left_data=[BaselineEMG_becips{1,trial_num}(:,1);FullEMG_becips{1,trial_num}(:,1)];
right_data=[BaselineEMG_becips{1,trial_num}(:,2);FullEMG_becips{1,trial_num}(:,2)];
% 绘制左手和右手数据在同一图中
plot( left_data, 'r', 'DisplayName', 'Left Hand');
hold on;
plot( right_data, 'b', 'DisplayName', 'Right Hand');
hold on;

% 设置图表标题和标签
title('Beicp');
xlabel('Time (s)');
ylabel('Amplitude');
% 添加图例
legend;

hold off;

%%

%test FDI
left_data=[BaselineEMG{1,trial_num}(:,1);FullEMG{1,trial_num}(:,1)];
right_data=[BaselineEMG{1,trial_num}(:,2);FullEMG{1,trial_num}(:,2)];
% 绘制左手和右手数据在同一图中
plot( left_data, 'r', 'DisplayName', 'Left Hand');
hold on;
plot( right_data, 'b', 'DisplayName', 'Right Hand');
hold on;

% 设置图表标题和标签
title('FDI');
xlabel('Time (s)');
ylabel('Amplitude');
% 添加图例
legend;

hold off;