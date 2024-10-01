%%  6.1 EMG analysis
% Author: Liang Tong
% Date: 12/7/2024 

%to update list
%1.function all the code to make main script simplier
%% Toolbox requirements: 
clc
clear all
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2'));% current folder

%% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);


%% 1. Load rl_BB pre-cue-188
% col 1: subject num
% col 2: block num
% col 3: contrast (high or low contrast)
% col 4: muscle (1=FDI, 2=BCP)
% col 5: trial outcome 345 is useless(1=correct, 2=error, 4=no response, 3 too early, 6 slow, 5 wrong muscle)
% col 6: RT in sec
% col 7: evshowtime =[]   % 0.2sec additional time after response
% col 8: participant response, 1 = left, 2 = right 0= on response
% col 9: correct response, 1 = left, 2 = right
% col 10: first tilt used for SSVEP

%new EMG onsite data created
%col 11 是更新后的 EMG onsite time
%col 12 是与原来时间的提前量
%col 13 是左手brust计数
%col 14 是右手brust计数
%col 15 是是否更新
addpath(genpath(exp.finalpath));

% load TL_ALL_include_EMG_updated.mat% EMG missing BCP totEMG
% 
% load TL_ALL_include_EMG% in EMGfolder

load([exp.behpath exp.name 'EMG_brusts_update'])

load([exp.behpath exp.name 'EMG_Vaild_check'])%trial目检有效性
load([exp.behpath exp.name 'EMGonsite_Lock_time'])%更新后时间
%% ，只要12 和valid
%1. 提前量对比
idx_R_FDI=find(EMGonsite_Valid==1 &AllBehaviour_EMGonsite(:,4)==4 &  (AllBehaviour_EMGonsite(:,5)== 1|   AllBehaviour_EMGonsite(:,5)== 2) );
idx_R_BCP=find(EMGonsite_Valid==1 &AllBehaviour_EMGonsite(:,4)==44 &  (AllBehaviour_EMGonsite(:,5)== 1|   AllBehaviour_EMGonsite(:,5)== 2) );

data1=AllBehaviour_EMGonsite(idx_R_FDI,6)-EMGonsite_lock_time(idx_R_FDI);%看提前量
data2=AllBehaviour_EMGonsite(idx_R_BCP,6)-EMGonsite_lock_time(idx_R_BCP);
%%
%2. EMG onsite 图 FDI vsBCP 更新后EMGonsite
idx_R_FDI=find(EMGonsite_Valid==1 &AllBehaviour_EMGonsite(:,4)==4 &  (AllBehaviour_EMGonsite(:,5)== 1|   AllBehaviour_EMGonsite(:,5)== 2) );
idx_R_BCP=find(EMGonsite_Valid==1&AllBehaviour_EMGonsite(:,4)==44 &  (AllBehaviour_EMGonsite(:,5)== 1|   AllBehaviour_EMGonsite(:,5)== 2) );


data1=EMGonsite_lock_time(idx_R_FDI);
data2=EMGonsite_lock_time(idx_R_BCP);
%%
%3. 数据质量对比 动了左手但右边的partial brust的个数
idx_R_FDI=find(AllBehaviour_EMGonsite(:,8)==2 &EMGonsite_Valid==1&AllBehaviour_EMGonsite(:,13)<8 &AllBehaviour_EMGonsite(:,14)<8 &AllBehaviour_EMGonsite(:,4)==4 &  (AllBehaviour_EMGonsite(:,5)== 1|   AllBehaviour_EMGonsite(:,5)== 2) );
idx_R_BCP=find(AllBehaviour_EMGonsite(:,8)==2 &EMGonsite_Valid==1&AllBehaviour_EMGonsite(:,13)<8 &AllBehaviour_EMGonsite(:,14)<8&AllBehaviour_EMGonsite(:,4)==44 &  (AllBehaviour_EMGonsite(:,5)== 1|   AllBehaviour_EMGonsite(:,5)== 2) );

data1=AllBehaviour_EMGonsite(idx_R_FDI,13);%看提前量
data2=AllBehaviour_EMGonsite(idx_R_BCP,13);

%%
% Create a new figure
figure;
nbins=30;
% Plot histogram for data1
h1 = histogram(data1,nbins, 'Normalization', 'probability');%'cdf'); % Normalized to show probability
hold on; % Retain current plot when adding new plots

% Plot histogram for data2
h2 = histogram(data2,nbins, 'Normalization', 'probability');%'probability');'cdf'); %

% Customize appearance
h1.FaceAlpha = 0.5; % Set transparency for the first histogram
h2.FaceAlpha = 0.5; % Set transparency for the second histogramh1
% h1.BinWidth = 1;    % Adjust bin width (optional)
% h2.BinWidth = 1;    % Adjust bin width (optional)

% Add legend and labels
legend('FDI', 'BCP');
xlabel('Value');
ylabel('Probability');

% Hold off to stop adding to the plot
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary statistics for data1
mean_data1 = mean(data1);
median_data1 = median(data1);
std_data1 = std(data1);
min_data1 = min(data1);
max_data1 = max(data1);
iqr_data1 = iqr(data1);

% Summary statistics for data2
mean_data2 = mean(data2);
median_data2 = median(data2);
std_data2 = std(data2);
min_data2 = min(data2);
max_data2 = max(data2);
iqr_data2 = iqr(data2);

% Display the summary statistics
disp(' ');
disp('Summary Statistics for Data1:');
disp(['Mean: ', num2str(mean_data1)]);
disp(['Median: ', num2str(median_data1)]);
disp(['Standard Deviation: ', num2str(std_data1)]);
disp(['Min: ', num2str(min_data1)]);
disp(['Max: ', num2str(max_data1)]);
disp(['Interquartile Range (IQR): ', num2str(iqr_data1)]);

disp(' ');

disp('Summary Statistics for Data2:');
disp(['Mean: ', num2str(mean_data2)]);
disp(['Median: ', num2str(median_data2)]);
disp(['Standard Deviation: ', num2str(std_data2)]);
disp(['Min: ', num2str(min_data2)]);
disp(['Max: ', num2str(max_data2)]);
disp(['Interquartile Range (IQR): ', num2str(iqr_data2)]);
% 统计学t检验
disp(' ');
% Perform two-sample t-test (unpaired, assuming unequal variance by default)
[h, p, ci, stats] = ttest2(data1, data2);

% Display results
%disp(['Hypothesis Test Result (h): ', num2str(h)]);
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: ', num2str(ci(1)), ' to ',num2str(ci(2))]);
%disp(['Test Statistic: ', num2str(stats.tstat)]);
%disp(['Degrees of Freedom: ', num2str(stats.df)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画箱型图看有没有奇怪分布
% Create a combined dataset for boxplotting
data_combined = [data1(:); data2(:)];

% Create a grouping variable to distinguish between the two datasets
group = [ones(size(data1(:))); 2 * ones(size(data2(:)))];

% Create a new figure
figure;

% Plot the boxplot
boxplot(data_combined, group, 'Labels', {'FDI', 'BCP'});

hold on; % Retain current plot for adding scatter points

% Add jittered scatter points for data1
x1 = ones(size(data1)) + (rand(size(data1)) - 0.5) * 0.2; % Jitter along the x-axis
scatter(x1, data1, 'filled', 'MarkerFaceAlpha', 0.05); % Plot jittered points for data1

% Add jittered scatter points for data2
x2 = 2 * ones(size(data2)) + (rand(size(data2)) - 0.5) * 0.2; % Jitter along the x-axis
scatter(x2, data2, 'filled', 'MarkerFaceAlpha', 0.05); % Plot jittered points for data2

% Customize appearance
ylabel('Values');
title('Boxplot with Jittered Scatter Points for Data1 and Data2');

hold off; % Stop adding to the plot