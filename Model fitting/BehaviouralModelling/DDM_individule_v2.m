%%--------------------------------------------------------------------------------------------
%%  1. DDM fitting script for individual calculation
% Author: Liang Tong
% Date: 12/3/2025


%%
clc
clear
minRT = 0.1;% can't be faster
maxRT = 2;
cohlevels = [0.07 0.14 0.07 0.14];
% read the raw behavioural data in:
data =readtable(['TL_ALL_ForR_1.csv']);
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
    if data_arr(i,3)==0.07 && data_arr(i,4)==4 % 4 is FDI
        temp_condition(i)=1;%1 = FDI_low
    elseif data_arr(i,3)==0.14 && data_arr(i,4)==4
        temp_condition(i)=2;%1 = FDI_high
    elseif data_arr(i,3)==0.07 && data_arr(i,4)==44
        temp_condition(i)=3;%1 = BCP_low
    elseif data_arr(i,3)==0.14 && data_arr(i,4)==44
        temp_condition(i)=4;%1 = BCP_high
    end

    %筛选结果对错和miss col 2
    if data_arr(i,8)==data_arr(i,9) && data_arr(i,5)~= 4 && data_arr(i,5)~= 3 && data_arr(i,5)~= 5
        temp_outcome(i)=1;%correct
    elseif data_arr(i,8)~=data_arr(i,9) && data_arr(i,5)~= 4 && data_arr(i,5)~= 3 && data_arr(i,5)~= 5
        temp_outcome(i)=0;%wrong
    elseif data_arr(i,5)== 4 ||data_arr(i,5)== 3||data_arr(i,5)== 5
        temp_outcome(i)=2;%miss
    end



end

data_allsubj(:,1)= temp_condition; %col 1: condition
data_allsubj(:,2)= temp_outcome; % col 2: trial outcome (1=correct, 0=error, 2=miss if relevant)
data_allsubj(:,3)= data_arr(:,6);% col 3: RT in sec (check it's not ms or EEG samples)
data_allsubj(:,4)=data_arr(:,1);% sub number

figure (1); hist(data_allsubj(:,3),[0:.001:5]); xlim([-0.1 4.5]); xlabel('RT') % note 1-ms resolution

miss = find(data_allsubj(:,3)>maxRT);
data_allsubj(miss,2)=maxRT;
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
datsum = makeSummaryStructure(data_allsubj,qps,4); % the 4 indicates the column with a subject indicator, 

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
    datsum = makeSummaryStructure(sub_data, qps, 0); % Compute summary per subject
    
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

% Define parameter ranges: 
par_range = [0 0.5;... % tnd 
            0 0.5;... % d(1) drift rate for lowest coherence FDI
            0 0.5;... % d(2) drift rate for high coherence FDI
            0 0.5;... % d(3) drift rate for lowest coherence BCP
            0 0.5;... % d(4)  drift rate for high coherence BCP
            0 0.5;... % b bound
            -0.25 0.5]; % y bound collapse slope - allowing it to increase over time, if that's what it needs to do to fit data



modelfn = 'simul_2AFCddm'; % name of function that will be called by G2 to generate predicted data
modelnames = {'ConstBnd_1d','ConstBnd_4d','CollBnd_1d','CollBnd_4d'};

% for each, designate which will be free to vary to fit data
Selvecs = {[1 1 0 1 0 1 0], [1 1 1 1 1 1 0], [1 1 0 1 0 1 1], [1 1 1 1 1 1 1]};
pm_if_fixed = [nan nan nan nan nan nan 0]; % order: [tnd d(x4) b y] any parameter that needs to be fixed at a certain value if not free to vary? 
N=1000; % the more trials simulated per condition, the more accurate the G^2, AIC etc, and less specific to the happenstance of the particular noise generating by the rng. Keeping it relatively low for this tutorial so doesn't take all day




tic

unique_subs = unique(data_allsubj(:,4));
n_subs = length(unique_subs);

results = struct(); % 初始化存储结果的结构体

for s = 1:n_subs
    sub_id = unique_subs(s);
    sub_data = data_allsubj(data_allsubj(:,4) == sub_id, :);
    datsum = makeSummaryStructure(sub_data, qps, 0); % 计算当前被试的数据摘要

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
            bestpm(m,5) = bestpm(m,2) * cohlevels(4) / cohlevels(3);
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

save('DDM_Fitting_individulexp1_Results.mat', 'results'); % 保存所有被试的拟合结果

% %% 初始化存储 AIC 和 BIC 的数组
% sub_ids = [];
% model_names = {};
% AIC_values = [];
% BIC_values = [];
% 
% % 遍历每个被试的结果
% for s = 1:length(results)
%     sub_id = results(s).sub_id; % 提取被试 ID
%     for m = 1:length(results(s).model)
%         model_name = results(s).model(m).name; % 提取模型名称
%         AIC = results(s).model(m).AIC; % 提取 AIC
%         BIC = results(s).model(m).BIC; % 提取 BIC
% 
%         % 存储数据
%         sub_ids = [sub_ids; sub_id];
%         model_names = [model_names; model_name];
%         AIC_values = [AIC_values; AIC];
%         BIC_values = [BIC_values; BIC];
%     end
% end
% 
% % 创建表格
% AIC_BIC_Table = table(sub_ids, model_names, AIC_values, BIC_values, ...
%     'VariableNames', {'Subject_ID', 'Model', 'AIC', 'BIC'});
% 
% % 显示结果
% disp(AIC_BIC_Table);

%%
% 初始化存储最优模型的数组
% best_models = cell(length(results), 2); % 第一列: Subject ID, 第二列: Best Model Name
% 
% % 遍历每个被试
% for s = 1:length(results)
%     sub_id = results(s).sub_id; % 获取被试 ID
%     min_AIC = inf; % 初始化最小 AIC
%     best_model = ''; % 初始化最优模型名称
% 
%     % 遍历该被试的所有模型
%     for m = 1:length(results(s).model)
%         current_AIC = results(s).model(m).AIC; % 获取 AIC 值
%         if current_AIC < min_AIC % 发现更小的 AIC
%             min_AIC = current_AIC;
%             best_model = results(s).model(m).name; % 更新最优模型名称
%         end
%     end
% 
%     % 存储最优模型信息
%     best_models{s, 1} = sub_id;
%     best_models{s, 2} = best_model;
% end
% 
% % 创建表格
% Best_Model_Table = cell2table(best_models, 'VariableNames', {'Subject_ID', 'Best_Model'});
% 
% % 显示结果
% disp(Best_Model_Table);

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
    datsum = makeSummaryStructure(sub_data, qps, 0); % 计算当前被试的数据摘要
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
    simdatsum = makeSummaryStructure(simdata, qps, 0);
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

%% 
 %计算所有被试的simdatsum和datsum的平均值
avg_simdatsum.n = zeros(size(simdatsum_all{1}.n));
avg_simdatsum.q = zeros(size(simdatsum_all{1}.q));
avg_simdatsum.qp = zeros(size(simdatsum_all{1}.qp));

avg_datsum.n = zeros(size(datsum_all{1}.n));
avg_datsum.q = zeros(size(datsum_all{1}.q));
avg_datsum.qp = zeros(size(datsum_all{1}.qp));

for s = 1:n_subs
    avg_simdatsum.n = avg_simdatsum.n + simdatsum_all{s}.n;
    avg_simdatsum.q = avg_simdatsum.q + simdatsum_all{s}.q;
    avg_simdatsum.qp = avg_simdatsum.qp + simdatsum_all{s}.qp;
    
    avg_datsum.n = avg_datsum.n + datsum_all{s}.n;
    avg_datsum.q = avg_datsum.q + datsum_all{s}.q;
    avg_datsum.qp = avg_datsum.qp + datsum_all{s}.qp;
end

avg_simdatsum.n = avg_simdatsum.n / n_subs;
avg_simdatsum.q = avg_simdatsum.q / n_subs;
avg_simdatsum.qp = avg_simdatsum.qp / n_subs;

avg_datsum.n = avg_datsum.n / n_subs;
avg_datsum.q = avg_datsum.q / n_subs;
avg_datsum.qp = avg_datsum.qp / n_subs;


  figure
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
    % title(['Subject ' num2str(sub_id) ', Model: ' best_model_name]);
    legend({'FDI_low','FL', 'FDI_high', 'FH', 'BCP_low','BL', 'BCP_high','BH'});
    hold off;
%% 总模型训练
datsum = makeSummaryStructure(data_allsubj,qps,4);

modelfn = 'simul_2AFCddm'; % name of function that will be called by G2 to generate predicted data
modelnames = {'ConstBnd_1d','ConstBnd_4d','CollBnd_1d','CollBnd_4d'};
% for each, designate which will be free to vary to fit data
Selvecs = {[1 1 0 1 0 1 0], [1 1 1 1 1 1 0], [1 1 0 1 0 1 1], [1 1 1 1 1 1 1]};
pm_if_fixed = [nan nan nan nan nan nan 0]; % order: [tnd d(x4) b y] any parameter that needs to be fixed at a certain value if not free to vary? Just y=0 for const bound
N=1000; % the more trials simulated per condition, the more accurate the G^2, AIC etc, and less specific to the happenstance of the particular noise generating by the rng. Keeping it relatively low for this tutorial so doesn't take all day
tic
for m=1:length(modelnames)
    disp(modelnames{m})
    Sel = Selvecs{m};
    options = optimoptions('particleswarm','UseParallel',true, 'SwarmSize', 10*sum(Sel));
    [p,g2(m)] = particleswarm(@(p) G2(p,pm_if_fixed(Sel==0),Sel,datsum,qps,modelfn,N,seed), sum(Sel), par_range((Sel==1),1)',par_range((Sel==1),2)', options);
    bestpm(m,Sel==1) = p; bestpm(m,Sel==0) = pm_if_fixed(Sel==0); 
    if all(Sel(2:5)==[1 0 1 0])
            % In the case of two free drift rate parameter, have to fill in the other drift rates:
            bestpm(m,3) = bestpm(m,2) * cohlevels(2) / cohlevels(1);
            %bestpm(m,4) = bestpm(m,2) * cohlevels(3) / cohlevels(1);
            bestpm(m,5) = bestpm(m,2) * cohlevels(4) / cohlevels(3);
        end
    disp(['best parameters: ' num2str(bestpm(m,:))])
    numfreepm=sum(Sel);
    
    % Do a fresh simulation of many trials with new noise, for computing AIC/BIC:
    g2new(m) = G2(p,pm_if_fixed(Sel==0),Sel,datsum,qps,modelfn,20000,sum(clock)); % note fresh seed and 20,000 trials per condn
    AIC(m) = g2new(m) + 2*numfreepm;
    BIC(m) = g2new(m) + numfreepm*log(sum(datsum.n(:)));
end

Sel = [1 1 1 1 1 1 1]; % the 2rd model has lowest AIC and BIC


simdata = simul_2AFCddm(bestpm(4,Sel==1),pm_if_fixed(Sel==0),Sel,10000,seed);
simdatsum = makeSummaryStructure(simdata,qps,0);
  figure
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
    
    legend({'FDI_low','FL', 'FDI_high', 'FH', 'BCP_low','BL', 'BCP_high','BH'});
    hold off;

    %% 初始化 response probability 变量


% Define RT quantile levels
quantiles = [0.1 0.3 0.5 0.6 0.9];

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

        % Compute response probabilities
        total_responses = sum(sub_cond_data(:,2) ~= 2);  % Exclude missed trials
        correct_responses = sum(sub_cond_data(:,2) == 1);
        error_responses = sum(sub_cond_data(:,2) == 0);

        if total_responses > 0
            p_correct = correct_responses / total_responses;
            p_error = error_responses / total_responses;
        else
            p_correct = NaN;
            p_error = NaN;
        end

              % Compute quantiles for correct responses
        correct_RT = sub_cond_data(sub_cond_data(:,2) == 1, 3); % RTs of correct trials
        if ~isempty(correct_RT)
            correct_quantiles = quantile(correct_RT, quantiles);
            % correct_hist = histcounts(correct_RT, [min(correct_RT)-eps, correct_quantiles, max(correct_RT)+eps]);
            % correctE_hist = histcounts(sub_cond_data(sub_cond_data(:,2) ~= 2, 3), [min(sub_cond_data(sub_cond_data(:,2) ~= 2, 3))-eps, correct_quantiles, max(sub_cond_data(sub_cond_data(:,2) ~= 2, 3))+eps]);
            % 
            % correct_prob_intervals = correct_hist ./ correctE_hist;
        else
            correct_quantiles = NaN(size(quantiles));
            % correct_prob_intervals = NaN(1, length(quantiles)+1);
        end

        % Compute quantiles for error responses
        error_RT = sub_cond_data(sub_cond_data(:,2) == 0, 3); % RTs of error trials
        if ~isempty(error_RT)
            error_quantiles = quantile(error_RT, quantiles);
            % error_hist = histcounts(error_RT, [min(error_RT)-eps, error_quantiles, max(error_RT)+eps]);
            % error_prob_intervals = error_hist / sum(error_hist);
        else
            error_quantiles = NaN(size(quantiles));
            % error_prob_intervals = NaN(1, length(quantiles)+1);
        end

        % Store results in a table
        result_table = [result_table; sub_id, condition, p_correct, p_error, ...
                        correct_quantiles, error_quantiles];
    end
end

% Convert to table format for better visualization
column_names = {'Subject', 'Condition', 'P_Correct', 'P_Error', ...
                'Q10_Correct', 'Q30_Correct', 'Q50_Correct', 'Q60_Correct', 'Q90_Correct', ...
                'Q10_Error', 'Q30_Error', 'Q50_Error', 'Q60_Error', 'Q90_Error'};

result_tables = array2table(result_table, 'VariableNames', column_names);



result_table(1,3) 
result_table(1,5:9)


figure
hold on
for c=1:4
plot(result_table(c,3),result_table(c,5:9),'o', 'Color', col(c,:)) % correct response probility
plot(result_table(c,4),result_table(c,10:14),'o', 'Color', col(c,:)) % error response probility

end

  xlabel('Response Probbility');
    ylabel('Reaction Time Quanitile(s)');
    
    legend({'FDI_low', 'FDI_high', 'BCP_low', 'BCP_high'});
hold off