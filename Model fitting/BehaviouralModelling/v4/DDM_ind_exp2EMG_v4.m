%%--------------------------------------------------------------------------------------------
%%  1. DDM fitting script for individual calculation
% Author: Liang Tong
% Date: 24/3/2025


%%
clc
clear
minRT = 0.1;% can't be faster
maxRT = 2.5;
cohlevels = [0.05 0.1 0.05 0.1];
% read the raw behavioural data in:
data =readtable(['TLb_ALL_ForR_exp2.csv']);
%example=load('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis\Model fitting\Modeling_G2_exDots\data_allsubj.mat');
% col 1: subject num
% col 2: block num
% col 3: contrast (high or low contrast)
% col 4: muscle (4=FDI, 44=BCP)
% col 5: trial outcome (1=correct, 2=error, 4=no response, 3 too early, 6 slow, 5 wrong muscle)
% col 6: RT in sec
% col 7: evshowtime =[]   % 0.2sec additional time after response
% col 8: participant response, 1 = left, 2 = right 0= on response
% col 9: correct response, 1 = left, 2 = right
% col 10: first tilt used for SSVEP
data_arr= table2array(data);
%%
% change to required format
% the behavioural data matrix should contain one row for each trial, and at least these 3 columns:
% col 1: condition (integers referring to different conditions. （1 = FDI_low，
% 2=FDI_high 3=BCP_Low 4=BCP_high)
% col 2: trial outcome (1=correct, 0=error, 2=miss if relevant)
% col 3: RT in sec (check it's not ms or EEG samples)

temp_condition=999*ones(size(data_arr,1),1);
temp_outcome=999*ones(size(data_arr,1),1);

for i=1: size(data_arr,1)

    %筛选条件 col 1
    if data_arr(i,3)==0.05 && data_arr(i,4)==4 % 4 is FDI
        temp_condition(i)=1;%1 = FDI_low
    elseif data_arr(i,3)==0.1 && data_arr(i,4)==4
        temp_condition(i)=2;%1 = FDI_high
    elseif data_arr(i,3)==0.05 && data_arr(i,4)==44
        temp_condition(i)=3;%1 = BCP_low
    elseif data_arr(i,3)==0.1 && data_arr(i,4)==44
        temp_condition(i)=4;%1 = BCP_high
    else
        temp_condition(i)=nan;
    end

    %筛选结果对错和miss col 2
    if data_arr(i,8)==data_arr(i,9) && data_arr(i,5)~= 4 && data_arr(i,5)~= 3 && data_arr(i,5)~= 5
        temp_outcome(i)=1;%correct
    elseif data_arr(i,8)~=data_arr(i,9) && data_arr(i,5)~= 4 && data_arr(i,5)~= 3 && data_arr(i,5)~= 5
        temp_outcome(i)=0;%wrong
    elseif data_arr(i,5)== 4 
        temp_outcome(i)=2;%miss
    elseif data_arr(i,5)== 3||data_arr(i,5)== 5
        temp_outcome(i)=nan;%delete
    end



end

data_allsubj(:,1)= temp_condition; %col 1: condition
data_allsubj(:,2)= temp_outcome; % col 2: trial outcome (1=correct, 0=error, 2=miss if relevant)
data_allsubj(:,3)= data_arr(:,11);% col 3: RT in sec (check it's not ms or EEG samples)
data_allsubj(:,4)=data_arr(:,1);% sub number

figure (1); hist(data_allsubj(:,3),[0:.001:5]); xlim([-0.1 4.5]); xlabel('RT') % note 1-ms resolution

miss = find(data_allsubj(:,3)>maxRT);
data_allsubj(miss,2)=2;
data_allsubj(miss,3)=nan;
%figure (2); hist(data_allsubj(:,3),[0:.001:5]); xlim([-0.1 2]); xlabel('RT') % note 1-ms resolution
miss = find(data_allsubj(:,3)<minRT);
data_allsubj(miss,2)=nan;
data_allsubj(miss,3)=nan;
figure (2); hist(data_allsubj(:,3),[0:.001:5]); xlim([-0.1 2.5]); xlabel('RT') % note 1-ms resolution
%% Plot RT histogram for each individual in a 2x3 subplot
unique_subs = unique(data_allsubj(:,4));
num_subs = length(unique_subs);
num_rows = 2;
num_cols = 3;
figure;
for s = 1:num_subs
    sub_id = unique_subs(s);
    sub_idx = find(data_allsubj(:,4) == sub_id);
    subplot(num_rows, num_cols, s);
    hist(data_allsubj(sub_idx,3), [0:.001:5]);
    xlim([-0.1 3]);
    xlabel('RT');
    ylabel('Frequency');
    title(['Subject ', num2str(sub_id)]);
end
%%
% Make the data summary:
% First define the RT quantiles from which we are going to compute the trial proportions for real and simulated data:
qps = [.1 .3 .5 .7 .9]; % these are quite standard
datsum = Copy_of_makeSummaryStructure(data_allsubj,qps,4); % the 4 indicates the column with a subject indicator, 

plotdatsum(datsum) % th
% 
% %% Compute summary data per subject
% qps = [.1 .3 .5 .7 .9]; % Define standard RT quantiles
% unique_subs = unique(data_allsubj(:,4));
% n_subs = length(unique_subs);
% 
% figure;
% for s = 1:n_subs
%     sub_id = unique_subs(s);
%     sub_data = data_allsubj(data_allsubj(:,4) == sub_id, :);
%     datsum = makeSummaryStructure(sub_data, qps, 0); % Compute summary per subject
% 
%     subplot(2, 3, s);
%     hold on;
%     col = lines(4);
%     lw = [1 3];
%     for ce = 1:2 % error and correct
%         for c = 1:4
%             plot((datsum.q(1:end-1,ce,c)+datsum.q(2:end,ce,c))/2, cumsum(datsum.qp(:,ce,c)),'.-', 'Color', col(c,:), 'LineWidth', lw(ce));
%         end
%     end
%     cohlevels = ["FDI_low", "FDI_high", "BCP_low", "BCP_high"];
%     legend(cohlevels);
%     title(['Subject ', num2str(sub_id)]);
%     hold off;
% end
%% Compute summary data per subject
qps = [.1 .3 .5 .7 .9]; % Define standard RT quantiles
unique_subs = unique(data_allsubj(:,4));
n_subs = length(unique_subs);

figure;
for s = 1:n_subs
    sub_id = unique_subs(s);
    sub_data = data_allsubj(data_allsubj(:,4) == sub_id, :);
    datsum = Copy_of_makeSummaryStructure(sub_data, qps, 0); % Compute summary per subject
    
    subplot(2, 3, s);
    hold on;
    col = lines(4);
    lw = [1 3];
    for ce = 1:2 % error and correct
        for c = 1:4
            plot((datsum.q(1:end-1,ce,c)+datsum.q(2:end,ce,c))/2, datsum.qp(:,ce,c),'.-', 'Color', col(c,:), 'LineWidth', lw(ce));
        end
    end
    cohlevels_name = ["FDI_low", "FDI_high", "BCP_low", "BCP_high"];
    legend(cohlevels_name);
    title(['Subject ', num2str(sub_id)]);
    hold off;
end
%% DDM modeling
seed = sum(clock); % pick one rng seed and use it for ALL simulations to minimise the jittering about
% % let's just set one drift rate for the lowest coherence and let the function scale up the rest
% Sel = [1 1 0 0 0 1 1];  
% pm = [0.2 0.05 0 0 0 0.08 0]; % order: [tnd d(x4) b y]
% simdata = simul_2AFCddm(pm(Sel==1),pm(Sel==0),Sel,1000,seed);
% simdatsum = makeSummaryStructure(simdata,qps,0);
% plotdatsum(simdatsum)
% 
% % The full parameter vector pm for this model is as follows:
% % tnd = nondecision time
% d(1) = drift rate for coherence level 1
% d(2) = drift rate for coherence level 2
% d(3) = drift rate for coherence level 3
% d(4) = drift rate for coherence level 4
% b = boundary height (separation = twice that) at time=0. UPPER BOUND IS CORRECT (convention)
% y = slope of bound collapse in units of /sec
% y2 for becips

% Define parameter ranges: 
par_range = [0 0.5;... % tnd 
            0 0.5;... % d(1) drift rate for lowest coherence FDI
            0 0.5;... % d(2) drift rate for high coherence FDI
            0 0.5;... % d(3) drift rate for lowest coherence BCP
            0 0.5;... % d(4)  drift rate for high coherence BCP
            0 0.5;... % b bound for FDI
            0 0.5;... % b bound for BCP
            -0.25 0.5;...% FDI y bound
            -0.25 0.5]; % y bound collapse slope - allowing it to increase over time, if that's what it needs to do to fit data



modelfn = 'Copy_of_simul_2AFCddm'; % name of function that will be called by G2 to generate predicted data
modelnames = {'ConstBnd_1d','ConstBnd_4d','CollBnd_1d','CollBnd_4d',...
                'FcollBnd_1d','FcollBnd_4d','BcollBnd_1d','BcollBnd_4d',...
                'ConstBnd_0d','CollBnd_0d','FcollBnd_0d','BcollBnd_0d'...
                'ConstBnd_0.5d','CollBnd_0.5d','FcollBnd_0.5d','BcollBnd_0.5d'...
                };

% for each, designate which will be free to vary to fit data
Selvecs = {[1 1 0 1 0 1 1 0 0], [1 1 1 1 1 1 1 0 0], [1 1 0 1 0 1 1 1 1], [1 1 1 1 1 1 1 1 1],...
    [1 1 0 1 0 1 1 1 0], [1 1 1 1 1 1 1 1 0], [1 1 0 1 0 1 1 0 1], [1 1 1 1 1 1 1 0 1],...
    [1 1 0 0 0 1 1 0 0],  [1 1 0 0 0 1 1 1 1], [1 1 0 0 0 1 1 1 0], [1 1 0 0 0 1 1 0 1]...
    [1 1 1 0 0 1 1 0 0],  [1 1 1 0 0 1 1 1 1], [1 1 1 0 0 1 1 1 0], [1 1 1 0 0 1 1 0 1]...
    };
pm_if_fixed = [nan nan nan nan nan nan nan 0 0]; % order: [tnd d(x4) b y] any parameter that needs to be fixed at a certain value if not free to vary? 
N=5000; % the more trials simulated per condition, the more accurate the G^2, AIC etc, and less specific to the happenstance of the particular noise generating by the rng. Keeping it relatively low for this tutorial so doesn't take all day




tic

unique_subs = unique(data_allsubj(:,4));
n_subs = length(unique_subs);

results = struct(); % 初始化存储结果的结构体

for s = 1:n_subs
    sub_id = unique_subs(s);
    sub_data = data_allsubj(data_allsubj(:,4) == sub_id, :);
    datsum = Copy_of_makeSummaryStructure(sub_data, qps, 0); % 计算当前被试的数据摘要

    for m=1:length(modelnames)
        disp(['Subject ' num2str(sub_id) ', Model: ' modelnames{m}])
        Sel = Selvecs{m};
        options = optimoptions('particleswarm','UseParallel',true, 'SwarmSize', 10*sum(Sel));
        
        [p, g2] = particleswarm(@(p) G2(p, pm_if_fixed(Sel==0), Sel, datsum, qps, modelfn, N, seed), sum(Sel), par_range((Sel==1),1)', par_range((Sel==1),2)', options);
        
        bestpm(m,Sel==1) = p; 
        bestpm(m,Sel==0) = pm_if_fixed(Sel==0); 
        
        if all(Sel(2:5)==[1 0 1 0])
            % In the case of two free drift rate parameter, have to fill in the other drift rates:
            bestpm(m,3) = bestpm(m,2) * cohlevels(2) / cohlevels(1);
            %bestpm(m,4) = bestpm(m,2) * cohlevels(3) / cohlevels(1);
            bestpm(m,5) = bestpm(m,4) * cohlevels(4) / cohlevels(3);

        elseif all(Sel(2:5)==[1 0 0 0])
           bestpm(m,3) = bestpm(m,2) * cohlevels(2) / cohlevels(1);
            bestpm(m,4) = bestpm(m,2) * cohlevels(3) / cohlevels(1);
            bestpm(m,5) = bestpm(m,2) * cohlevels(4) / cohlevels(1);

        elseif all(Sel(2:5)==[1 1 0 0])
           %bestpm(m,3) = bestpm(m,2) * cohlevels(2) / cohlevels(1);
            bestpm(m,4) = bestpm(m,2) ;
            bestpm(m,5) = bestpm(m,3) ;

        end
        
        disp(['Best parameters for Subject ' num2str(sub_id) ': ' num2str(bestpm(m,:))])
        numfreepm = sum(Sel);
        
        % Compute AIC/BIC with fresh simulations
        g2new = G2(p, pm_if_fixed(Sel==0), Sel, datsum, qps, modelfn, 20000, sum(clock));
        AIC = g2new + 2 * numfreepm;
        BIC = g2new + numfreepm * log(sum(datsum.n(:)));
        
        % 存储每个被试的数据
        results(s).sub_id = sub_id;
        results(s).model(m).name = modelnames{m};
        results(s).model(m).params = bestpm(m, :);
        results(s).model(m).G2 = g2;
        results(s).model(m).G2new = g2new;
        results(s).model(m).AIC = AIC;
        results(s).model(m).BIC = BIC;
    end
end

toc % 记录运行时间

save('DDM_Fitting_individulexp2EMG_Results_v2.mat', 'results'); % 保存所有被试的拟合结果
%%
%% Extract Data for Boxplots
n_models = length(modelnames);
n_subs = length(results);

AIC_values = NaN(n_subs, n_models);
BIC_values = NaN(n_subs, n_models);
G2_values = NaN(n_subs, n_models);

for s = 1:n_subs
    for m = 1:n_models
        if isfield(results(s).model(m), 'AIC')
            AIC_values(s, m) = results(s).model(m).AIC;
            BIC_values(s, m) = results(s).model(m).BIC;
            G2_values(s, m) = results(s).model(m).G2;
        end
    end
end

%% Boxplot - AIC, BIC, G²
figure;

% --- AIC ---
subplot(1,3,1);
hold on;
boxplot(AIC_values, 'Labels', modelnames, 'Symbol', 'r+'); % Boxplot with red outliers
%scatter(repmat(1:n_models, n_subs, 1) + (rand(n_subs, n_models) - 0.5) * 0.1, AIC_values(:), 30, 'k', 'filled'); % Jitter scatter
ylabel('AIC');
title('AIC Comparison');
grid on;

% --- BIC ---
subplot(1,3,2);
hold on;
boxplot(BIC_values, 'Labels', modelnames, 'Symbol', 'r+'); % Boxplot with red outliers
%scatter(repmat(1:n_models, n_subs, 1) + (rand(n_subs, n_models) - 0.5) * 0.1, BIC_values(:), 30, 'k', 'filled'); % Jitter scatter
ylabel('BIC');
title('BIC Comparison');
grid on;

% --- G² ---
subplot(1,3,3);
hold on;
boxplot(G2_values, 'Labels', modelnames, 'Symbol', 'r+'); % Boxplot with red outliers
%scatter(repmat(1:n_models, n_subs, 1) + (rand(n_subs, n_models) - 0.5) * 0.1, G2_values(:), 30, 'k', 'filled'); % Jitter scatter
ylabel('G²');
title('G² Comparison');
grid on;

hold off;


%% 显示平均参数
% Define Parameter Names
param_names = {'tnd', 'd1_FDI_low', 'd2_FDI_high', 'd3_BCP_low', 'd4_BCP_high', ...
               'b_FDI','b_BCP', 'y_FDI', 'y_BCP'}; % Update parameter names accordingly

% Extract Model Parameters and Compute Averages
n_models = length(modelnames);
n_subs = length(results);
num_params = length(results(1).model(1).params);

% Initialize parameter storage
param_matrix = NaN(n_subs, n_models, num_params);

for s = 1:n_subs
    for m = 1:n_models
        param_matrix(s, m, :) = results(s).model(m).params;
    end
end

% Compute the mean parameters across subjects
avg_params = squeeze(nanmean(param_matrix, 1)); % Average across subjects

% Convert results to table
avg_param_table = array2table(avg_params, 'VariableNames', param_names);
avg_param_table.Model = modelnames';
avg_param_table = movevars(avg_param_table, 'Model', 'Before', 'tnd');

% Display the averaged parameter table
disp(avg_param_table);

%%
%重画模型
figure;
tiledlayout(2,3); % 创建 2x3 子图布局

simdatsum_all = cell(n_subs,1); % 初始化存储所有被试的simdatsum
datsum_all = cell(n_subs,1); % 初始化存储所有被试的simdatsum

for s = 1:n_subs
    sub_id = unique_subs(s);
    disp(['Simulating for Subject ' num2str(sub_id)])

%原始数据
    sub_data = data_allsubj(data_allsubj(:,4) == sub_id, :);
    datsum = Copy_of_makeSummaryStructure(sub_data, qps, 0); % 计算当前被试的数据摘要
datsum_all{s}=datsum;
    
    % 获取该被试 AIC 最小的模型
    best_model_idx = -1;
    min_AIC = inf;
    best_model_name = '';
    
    for m = 1:length(results(s).model)
        if results(s).model(m).AIC < min_AIC
            min_AIC = results(s).model(m).AIC;
            best_model_idx = m;
            best_model_name = results(s).model(m).name;
        end
    end
    
    % 确保找到了最优模型
    if best_model_idx == -1
        disp(['No valid model found for Subject ' num2str(sub_id)])
        continue;
    end
    
    % 获取该被试最优模型的参数
    best_Sel = Selvecs{best_model_idx};
    best_params = results(s).model(best_model_idx).params;
    
    % 运行模拟
    simdata = simul_2AFCddm(best_params(best_Sel==1), pm_if_fixed(best_Sel==0), best_Sel, 10000, seed);
    simdatsum = Copy_of_makeSummaryStructure(simdata, qps, 0);
    simdatsum_all{s}=simdatsum;
    
    disp(['Computed simdatsum for Subject ' num2str(sub_id) ', Model: ' best_model_name]);
    
    % 画子图
    nexttile;
    hold on;
    col = lines(4);
    lw = [1 3]; % 线宽 - 错误薄，正确粗
    
    for ce = 1:2 % 错误和正确
        for c = 1:4
            plot((simdatsum.q(1:end-1,ce,c) + simdatsum.q(2:end,ce,c))/2, ...
                simdatsum.qp(:,ce,c), '.-', 'Color', col(c,:), 'LineWidth', lw(ce));
             plot((datsum.q(1:end-1,ce,c)+datsum.q(2:end,ce,c))/2, datsum.qp(:,ce,c),'o', 'Color', col(c,:), 'LineWidth', lw(ce));
        end
    end
    
    xlabel('Reaction Time (s)');
    ylabel('Probability Density');
    title(['Subject ' num2str(sub_id) ', Model: ' best_model_name]);
    legend({'FDI_low','FL', 'FDI_high', 'FH', 'BCP_low','BL', 'BCP_high','BH'});
    hold off;
end
%% 平均AIC最小
% Determine the Best Model (Lowest Average AIC)
n_models = length(modelnames);
n_subs = length(results);
AIC_values = NaN(n_subs, n_models);

for s = 1:n_subs
    for m = 1:n_models
        AIC_values(s, m) = results(s).model(m).AIC;
    end
end

% Compute average AIC per model
avg_AIC = nanmean(AIC_values, 1);
[~, best_model_idx] = min(avg_AIC);
best_model_name = modelnames{best_model_idx};

disp(['Selected model: ' best_model_name ' (Lowest average AIC)']);
%% Simulate & Plot Using the Best Model for All Subjects
figure;
tiledlayout(2,3); % 创建 2x3 子图布局

simdatsum_all = cell(n_subs,1); % 初始化存储所有被试的simdatsum
datsum_all = cell(n_subs,1); % 初始化存储所有被试的datsum

for s = 1:n_subs
    sub_id = unique_subs(s);
    disp(['Simulating for Subject ' num2str(sub_id) ' using model ' best_model_name]);

    % 原始数据
    sub_data = data_allsubj(data_allsubj(:,4) == sub_id, :);
    datsum = Copy_of_makeSummaryStructure(sub_data, qps, 0); % 计算当前被试的数据摘要
    datsum_all{s} = datsum;
    
    % 获取该被试最佳模型的参数
    best_Sel = Selvecs{best_model_idx};
    best_params = results(s).model(best_model_idx).params;
    
    % 运行模拟
    simdata = simul_2AFCddm(best_params(best_Sel==1), pm_if_fixed(best_Sel==0), best_Sel, 10000, seed);
    simdatsum = Copy_of_makeSummaryStructure(simdata, qps, 0);
    simdatsum_all{s} = simdatsum;
    
    disp(['Computed simdatsum for Subject ' num2str(sub_id)]);

    % 画子图
    nexttile;
    hold on;
    col = lines(4);
    lw = [1 3]; % 线宽 - 错误薄，正确粗
    
    for ce = 1:2 % 错误和正确
        for c = 1:4
            plot((simdatsum.q(1:end-1,ce,c) + simdatsum.q(2:end,ce,c))/2, ...
                simdatsum.qp(:,ce,c), '.-', 'Color', col(c,:), 'LineWidth', lw(ce));
            plot((datsum.q(1:end-1,ce,c) + datsum.q(2:end,ce,c))/2, ...
                datsum.qp(:,ce,c), 'o', 'Color', col(c,:), 'LineWidth', lw(ce));
        end
    end
    
    xlabel('Reaction Time (s)');
    ylabel('Probability Density');
    title(['Subject ' num2str(sub_id) ', Model: ' best_model_name]);
    legend({'FDI_low','FL', 'FDI_high', 'FH', 'BCP_low','BL', 'BCP_high','BH'});
    hold off;
end

% %%
% %显示所有参数
% % Extract Model Parameters and Compute Averages
% n_models = length(modelnames);
% n_subs = length(results);
% 
% % Initialize parameter storage
% param_matrix = NaN(n_subs, n_models, length(results(1).model(1).params));
% 
% for s = 1:n_subs
%     for m = 1:n_models
%         param_matrix(s, m, :) = results(s).model(m).params;
%     end
% end
% 
% % Compute the mean parameters across subjects
% avg_params = squeeze(nanmean(param_matrix, 1)); % Average across subjects
% 
% % Convert results to table
% avg_param_table = array2table(avg_params, 'VariableNames', strcat('Param_', string(1:size(avg_params,2))));
% avg_param_table.Model = modelnames';
% avg_param_table = movevars(avg_param_table, 'Model', 'Before', 'Param_1');
% 
% % Display the averaged parameter table
% disp(avg_param_table);

%% 初始化 response probability 变量

% Define RT quantile levels
quantiles = [0.1 0.3 0.5 0.7 0.9];

% Get unique subjects and conditions
unique_subs = unique(data_allsubj(:,4));  % Subject IDs
unique_conditions = unique(data_allsubj(:,1)); % Condition IDs

% Initialize result storage
result_table = [];

% Loop through each subject and condition
for s = 1:length(unique_subs)
    sub_id = unique_subs(s);
    for c = 1:length(unique_conditions)
        condition = unique_conditions(c);

        % Extract trials for this subject and condition
        sub_cond_data = data_allsubj(data_allsubj(:,4) == sub_id & data_allsubj(:,1) == condition, :);

        % Compute response counts
        %total_trials = size(sub_cond_data, 1);
        correct_responses = sum(sub_cond_data(:,2) == 1);
        error_responses   = sum(sub_cond_data(:,2) == 0);
        miss_responses    = sum(sub_cond_data(:,2) == 2);
        total_trials=correct_responses+error_responses+miss_responses; %%%%%这里改变总数定义

        % Compute response probabilities
        if total_trials > 0
            p_correct = correct_responses / total_trials;
            p_error   = error_responses / total_trials;
            p_miss    = miss_responses / total_trials;
        else
            p_correct = NaN;
            p_error   = NaN;
            p_miss    = NaN;
        end

        % Compute quantiles for correct responses
        correct_RT = sub_cond_data(sub_cond_data(:,2) == 1, 3);
        if ~isempty(correct_RT)
            correct_quantiles = quantile(correct_RT, quantiles);
        else
            correct_quantiles = NaN(size(quantiles));
        end

        % Compute quantiles for error responses
        error_RT = sub_cond_data(sub_cond_data(:,2) == 0, 3);
        if ~isempty(error_RT)
            error_quantiles = quantile(error_RT, quantiles);
        else
            error_quantiles = NaN(size(quantiles));
        end

        % Store results
        result_table = [result_table; sub_id, condition, p_correct, p_error, p_miss, ...
                        correct_quantiles, error_quantiles];
    end
end

% Convert to table
column_names = {'Subject', 'Condition', 'P_Correct', 'P_Error', 'P_Miss', ...
                'Q10_Correct', 'Q30_Correct', 'Q50_Correct', 'Q70_Correct', 'Q90_Correct', ...
                'Q10_Error', 'Q30_Error', 'Q50_Error', 'Q70_Error', 'Q90_Error'};

result_tables = array2table(result_table, 'VariableNames', column_names);



% 计算每个条件的 group-level 平均结果
group_avg_result_table = [];

for c = 1:length(unique_conditions)
    condition = unique_conditions(c);
    
    % 提取该条件下的所有行
    cond_rows = result_table(result_table(:,2) == condition, :);  % Condition列是第2列
    
    % 计算各列平均值（从 P_Correct 到最后一个 quantile）
    mean_vals = mean(cond_rows(:,3:end), 1, 'omitnan');  % Skip Subject, Condition
    
    % 组装：用 0 作为 group-level Subject ID
    group_avg_row = [0, condition, mean_vals];
    
    % 添加到 group 表格中
    group_avg_result_table = [group_avg_result_table; group_avg_row];
end

% 与原始 column_names 对齐
% column_names = {'Subject', 'Condition', 'P_Correct', 'P_Error', 'P_Miss', ...
%                 'Q10_Correct', 'Q30_Correct', 'Q50_Correct', 'Q70_Correct', 'Q90_Correct', ...
%                 'Q10_Error', 'Q30_Error', 'Q50_Error', 'Q70_Error', 'Q90_Error'};
% 
% group_avg_result_table = array2table(group_avg_result_table, 'VariableNames', column_names);

%% for all average
x= zeros(1,8);
y=zeros(5,8);

extract_table=group_avg_result_table;
x(1:4)= extract_table(1:4,4); %error1-4
x(5:8)= extract_table(1:4,3); %correct1-4

for i=1:5
    y(i,1:4)=extract_table(1:4,11+i-1); %Q10 error
    y(i,5:8)=extract_table(1:4,6+i-1); %Q10 correct
end
col = lines(4);
figure
hold on
col = lines(4);
for i=1:5
    plot(x(1:2),y(i,1:2), 'o-', 'Color', col(1,:))% fdi ERROR
    plot(x(5:6),y(i,5:6), 'o-', 'Color', col(1,:))
end
for i=1:5
    plot(x(3:4),y(i,3:4), '*-', 'Color', col(3,:))% bcp ERROR
    plot(x(7:8),y(i,7:8), '*-', 'Color', col(3,:))
end

xlabel('Response Probability');
ylabel('RT Quantile (s)');
%legend('Location','best');
title('Per-Condition RT Quantiles vs. Response Probability');
hold off
%% 单独画
% 提取所有被试编号
unique_subs = unique(result_table(:,1));
n_subs = length(unique_subs);

col = lines(4); % 保持颜色一致
quantile_labels = {'Q10', 'Q30', 'Q50', 'Q70', 'Q90'};

figure;
tiledlayout(2,3); % 2行3列的子图布局

for s = 1:n_subs
    sub_id = unique_subs(s);
    extract_table = result_table(result_table(:,1) == sub_id, :); % 提取当前被试数据
    
    % 初始化 x 和 y（x: response prob, y: quantile）
    x = zeros(1,8);
    y = zeros(5,8);
    
    % error: columns 4 = P_Error
    x(1:4) = extract_table(1:4,4); % 条件1-4的错误概率
    % correct: columns 3 = P_Correct
    x(5:8) = extract_table(1:4,3); % 条件1-4的正确概率
    
    for i = 1:5
        y(i,1:4) = extract_table(1:4,10+i); % Q10–Q90 Error: cols 11–15
        y(i,5:8) = extract_table(1:4,5+i);  % Q10–Q90 Correct: cols 6–10
    end

    nexttile; % 移动到下一个 subplot
    hold on

    for i = 1:5
        plot(x(1:2), y(i,1:2), 'o-', 'Color', col(1,:)) % FDI Error
        plot(x(5:6), y(i,5:6), 'o-', 'Color', col(1,:)) % FDI Correct
        plot(x(3:4), y(i,3:4), '*-', 'Color', col(2,:)) % BCP Error
        plot(x(7:8), y(i,7:8), '*-', 'Color', col(2,:)) % BCP Correct
    end

    xlabel('Response Probability');
    ylabel('RT Quantile (s)');
    title(['Subject ' num2str(sub_id)]);
    %ylim([0.1 1.5]); % 你可以根据需要调整范围
    hold off
end

legend({'FDI Error','FDI Correct','BCP Error','BCP Correct'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
%% miss 率
% 使用正确的表变量名
all_sub_ids = unique(result_tables.Subject);
all_conditions = unique(result_tables.Condition);

% 初始化遗漏概率矩阵（包括 group 平均，共 n+1 行）
miss_matrix = NaN(length(all_sub_ids)+1, 4); 

for s = 1:length(all_sub_ids)
    sub_id = all_sub_ids(s);
    for c = 1:length(all_conditions)
        cond = all_conditions(c);
        row_idx = result_tables.Subject == sub_id & result_tables.Condition == cond;
        if any(row_idx)
            miss_matrix(s, c) = result_tables.P_Miss(row_idx);
        end
    end
end

% 计算 group average（不含 NaN）
miss_matrix(end, :) = mean(miss_matrix(1:end-1,:), 'omitnan');

% 构建行名（含 Sub_x 和 Group_Avg）
row_names = [cellstr("Sub_" + string(all_sub_ids)); "Group_Avg"];

% 构建列名
var_names = {'FDI_Low', 'FDI_High', 'BCP_Low', 'BCP_High'};

% 转换为 table
miss_table = array2table(miss_matrix, 'VariableNames', var_names, 'RowNames', row_names);

% 显示结果
disp("Miss Probability per Subject and Group Average:");
disp(miss_table);

