function [EEG] = TLBF2_updateEMGtrigger(sub,exp)

% % Load converted data
EEG = pop_loadset([exp.filepath 'aac' exp.name num2str(sub) '.set']);

%找到证据开始点
index_trigger10x = find(ismember({EEG.event.type}, {'101', '102', '103', '104'}));
index_trigger10 = find(ismember({EEG.event.type}, '10'));

if length(index_trigger10x)~= 1024 && length(index_trigger10)~= 1024%检查数量
    error('length problem, please check');
    
end

%准备导入EMG数据
EMG_onsite=zeros(1024,1);
EMG_valid=zeros(1024,1);
load([exp.behpath exp.name 'EMGonsite_Lock_time'])%更新后时间
%根据sub写入时间
EMG_onsite=EMGonsite_lock_time((sub-1)*1024+1:(sub)*1024);
EMG_valid=EMGonsite_Valid((sub-1)*1024+1:(sub)*1024);
clear EMGonsite_lock_time;clear EMGonsite_Valid

%准备延后时间
fs = EEG.srate; % 获取采样率
latency_shift= (232+80) * (fs / 1000); % 结束前232秒一般为应答时间如无法更新


count=length(EEG.event);

for i = 1:1024

    if ~isnan(EMG_onsite(i)) %如有效则更新
    % 新建事件
    new_event.latency = EEG.event(index_trigger10x(i)).latency + EMG_onsite(i)* (fs ); % 在当前触发器之后的1个单位时间加1660
    new_event.type = '66';
    new_event.duration = 0;
    new_event.urevent=count+i;
    
    % 插入到EEG.event中
    EEG.event(end+1) = new_event;%EEG.event(end+1) = EEG.event(i);

    elseif isnan(EMG_onsite(i)) || EMG_valid(i)==0 %如无效则按照结束前232秒插入
    new_event.latency = EEG.event(index_trigger10(i)).latency - latency_shift; % 在当前触发器之后的1个单位时间加1660
    new_event.type = '66';
    new_event.duration = 0;
    new_event.urevent=count+i;
    
    % 插入到EEG.event中
    EEG.event(end+1) = new_event;
    else
        error('插入值问题，请检查');
        break

    end
%disp(i);
end

%重新排序
EEG = pop_editeventvals(EEG, 'sort', {'latency', [0]});

% %%保存和检查
% EEG = eeg_checkset( EEG );
% if exp.filter.lowerbound == 0.1
%     filename = ['f01' EEG.filename];
% elseif exp.filter.lowerbound == 0.01
%     filename = ['f' EEG.filename];
% end

% EEG.filename = filename;

%EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset

%filename = ['aac' exp.name num2str(sub)];
%EEG = pop_saveset(EEG, filename, exp.filepath);
% Display success message
disp(['EOL  added  successfully.']);
end