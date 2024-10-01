function [EEG1 EEG2] = TLBF2_baselineData(sub,exp,EEG, default)
% ver2: 27 Sep 2024
%加了一个EoL的lock

%% BASELINES:
% Calculate different baselines w
% PreStimulus: before lead-in dots --> reveals anticipatory build-up
% PrePulse: before first pulse --> should cancel anticipatory build-up
% Action-locked: around action (-5 +5)

% if contains(EEG.filename, 'SL') & ~contains(EEG.filename, 'RL')
%     EEG1 = EEG;
%     EEG2 = pop_loadset( [exp.filepath '/cICAtrial_riRL_' exp.filterLab 'atmc' (exp.name) '_P' num2str(sub) '.set'] );
% elseif contains(EEG.filename, 'RL')
%     EEG2 = EEG;
%     EEG1 = pop_loadset([exp.filepath '/cICAtrial_eyeB_rej_riP1_' exp.filterLab 'atmc' (exp.name) '_P' num2str(sub) '.set']);
% end

% %%Reject same trials that were rejected due to eyeblinks during pulse
% toReject = EEG1.eyeBlink_rej;
% EEG2.eyeBlink_rej = toReject;
% EEG2.sumData(toReject,:) = [];
%
% %Actually remove from data
% EEG2.data(:,:,EEG2.eyeBlink_rej) = [];
% EEG2.epoch(EEG2.eyeBlink_rej) = [];
% EEG2.trials = size(EEG2.data, 3);
%
% disp(['Percentage trials w/eyeblink during pulses: ' num2str(EEG1.eyeBlink_rej_perc)])
EEG1=EEG;%主要保存的b_数据
EEG2=EEG;%次要保存的bb_数据
EEG1_preStim = EEG1; %EEG2_preStim = EEG2;
EEG1_prePulse = EEG1; %EEG2_prePulse = EEG2;
EEG1_action = EEG1; %EEG2_action = EEG2;

%% PrePulse baseline (-0.18s before evidence 4 times SSVEP)
%这个为切割时0点前-188秒作为基线
%Calculate prePulse baseline - this should neutralise effects of
%anticipatory build-up

t.prePulse_baseline_idx = ([EEG1.times] >= -188 & [EEG1.times] <= 0);
t.prePulse_baseline = squeeze(mean(EEG1.data(:,t.prePulse_baseline_idx,:),2));

%% Action-locked baseline (all the trial)
% if contains(EEG.filename, 'RL') %单独找，用pre cue daseline
% t.actionBaseline_idx = ([EEG1.times] >= -5 & [EEG1.times] <= 5);
% t.actionBaseline = squeeze(mean(EEG1.data(:,t.actionBaseline_idx,:),2));
% end

%% Pre-stimulus baseline (-200 ms pre cue)
%一般都用cue'6'前-188ms作为baseline,但切割时不是'6'，所以要单独取值
% Identify 200ms preceding lead-in dots and use as preStim baseline

for e = 1:length(EEG1.epoch)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %pres cue的baseline
    %%%%小补丁
    % 假设输入是[EEG1.epoch(e).eventtype(1,:)]
    inputData = [EEG1.epoch(e).eventtype(1,:)];

    % 目标值
    if default(1)=='S'
        targetValue = '6';
    elseif default(1) =='R'
        targetValue = '6';
    elseif default(1) =='E'
        targetValue = '6';
    end

    % 查找目标值的位置并存储索引
    t.lidots_idx = [];
    for i = 1:length(inputData)
        if strcmp(inputData{i}, targetValue)
            t.lidots_idx = [t.lidots_idx, i];  % 如果找到目标值，将索引添加到idx数组
        end
    end

    % t.lidots_idx = cell2mat([EEG1.epoch(e).eventtype(1,:)]) == '6';
    t.lidots_time = cell2mat([EEG1.epoch(e).eventlatency(1,t.lidots_idx)]);
    t.preStim_baseline = [t.lidots_time-188, t.lidots_time]; %188is 4*SSVEP 21.5Hz
    t.preStim_baseline_idx = ([EEG1.times] >= t.preStim_baseline(1) & [EEG1.times] <= t.preStim_baseline(end));
    t.preStim_baseline = mean(EEG1.data(:,t.preStim_baseline_idx,e),2);

    EEG1.preStim_baseline(:,e) = t.preStim_baseline;
    EEG2.preStim_baseline(:,e) = t.preStim_baseline;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %pre cue的baseline 为RL单独找
    % t.actionBaseline_idx = ([EEG1.times] >= -5 & [EEG1.times] <= 5);
    % t.actionBaseline = squeeze(mean(EEG1.data(:,t.actionBaseline_idx,:),2));
    
    %%%%小补丁
    % 假设输入是[EEG1.epoch(e).eventtype(1,:)]
    inputData_RL = [EEG1.epoch(e).eventtype(1,:)];

    % 目标值

    targetValue_RL = '66'; %这里指的是EoL,本来粗略的用'10'前90ms估计（其实应该是232ms）;现统一跟新为'66'（EMG onsite time）


    % 查找目标值的位置并存储索引
    t.lidots_idx_RL = [];
    for i = 1:length(inputData_RL)
        if strcmp(inputData_RL{i}, targetValue_RL)
            t.lidots_idx_RL = [t.lidots_idx_RL, i];  % 如果找到目标值，将索引添加到idx数组
        end
    end

    % t.lidots_idx = cell2mat([EEG1.epoch(e).eventtype(1,:)]) == '6';
    t.lidots_time_RL = cell2mat([EEG1.epoch(e).eventlatency(1,t.lidots_idx_RL)]);
    t.actionBaseline = [t.lidots_time_RL-188, t.lidots_time_RL]; %188is 4*SSVEP 21.5Hz
    t.actionBaseline_idx = ([EEG1.times] >= t.actionBaseline(1) & [EEG1.times] <= t.actionBaseline(end));
    t.actionBaseline = mean(EEG1.data(:,t.actionBaseline_idx,e),2);

    EEG1.actionBaseline(:,e) = t.actionBaseline;
    EEG2.actionBaseline(:,e) = t.actionBaseline;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PreStim
    EEG1_preStim.data(:,:,e) = EEG1.data(:,:,e) - t.preStim_baseline; %Baseline data
    %EEG2_preStim.data(:,:,e) = EEG2.data(:,:,e) - t.preStim_baseline; %Baseline data

    %PrePulse
    EEG1_prePulse.data(:,:,e) = EEG1.data(:,:,e) - t.prePulse_baseline(:,e); %Baseline data
    %EEG2_prePulse.data(:,:,e) = EEG2.data(:,:,e) - t.prePulse_baseline(:,e); %Baseline data

    %     %Action
    EEG1_action.data(:,:,e) = EEG1.data(:,:,e) - t.actionBaseline; %Baseline data
    %EEG2_action.data(:,:,e) = EEG2.data(:,:,e) - t.actionBaseline(:,e); %Baseline data

end
EEG1.prePulse_baseline = t.prePulse_baseline;
%EEG1.preStim_baseline = t.preStimBaseline; 每个不同所以在loop里
EEG1.action_baseline = t.actionBaseline;
EEG2.prePulse_baseline = t.prePulse_baseline;
%EEG2.preStim_baseline = t.preStimBaseline; 每个不同所以在loop里
EEG2.action_baseline = t.actionBaseline;
% EEG2.prePulse_baseline = t.prePulse_baseline;
% EEG2.actionBaseline = t.actionBaseline;
%
EEG.nbchan
EEG1 = eeg_checkset( EEG1);
EEG2 = eeg_checkset( EEG2);

%By default, subtract:
if default(1)=='S'
    EEG1.data = EEG1_prePulse.data; %证据前0.2秒
    EEG2.data = EEG1_preStim.data;% cue前0.2秒
elseif default(1)=='R'
    EEG1.data =  EEG1_preStim.data;%结束前0.2秒%后期调整为EMG更准确的时间 
    EEG2.data =  EEG1_action.data;% cue前0.2秒

elseif default(1)=='E'
    EEG1.data =  EEG1_preStim.data;%结束前0.2秒%后期调整为EMG更准确的时间 
    EEG2.data =  EEG1_action.data;% cue前0.2秒
    % case 'action'
    %     EEG1.data = EEG1_action.data;
    %     EEG2.data = EEG2_action.data;
end

disp(size(EEG1.data));


EEG1.filename = ['b_' EEG1.filename];
EEG2.filename = ['bb_' EEG2.filename];


% EEG1.data(129:136,:,:) = []; %Drop last 8 channel;
% EEG2.data(129:136,:,:) = []; %Drop last 8 channel;

EEG1 = eeg_checkset( EEG1);
EEG2 = eeg_checkset( EEG2);
%
% if contains(EEG.filename, 'SL') & ~contains(EEG.filename, 'RL')
%     EEG = EEG1;
% elseif contains(EEG.filename, 'RL')
%     EEG = EEG2;
% end
%
EEG1= pop_saveset(EEG1, EEG1.filename, exp.filepath); %  for corrected
EEG2= pop_saveset(EEG2, EEG2.filename, exp.filepath); %  for corrected

% EEG2= pop_saveset(EEG2, ['b_' EEG2.filename], exp.filepath); %  for corrected


end