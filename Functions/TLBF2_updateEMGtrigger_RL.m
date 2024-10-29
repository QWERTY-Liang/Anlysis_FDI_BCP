function [EEG1] = TLBF2_updateEMGtrigger_RL(EEG,exp,sub)

% % Load converted data
%EEG = pop_loadset([exp.filepath 'aac' exp.name num2str(sub) '.set']);

%找到证据开始点
index_trigger10x = find(ismember({EEG.event.type}, {'101', '102', '103', '104'}));
index_trigger10 = find(ismember({EEG.event.type}, '10'));

if length(index_trigger10x)~= 1024 && length(index_trigger10)~= 1024%检查数量
    error('length problem, please check');
    
end


%准备导入EMG数据
EMG_onsite=zeros(1024,1);
EMG_valid=zeros(1024,1);
load([exp.behpath exp.name '_All_include_EMG'])%更新后时间
%根据sub写入时间
EMG_onsite=AllBehaviour((sub-1)*1024+1:(sub)*1024,6);
EMG_end=AllBehaviour((sub-1)*1024+1:(sub)*1024,7);% trial 结束时间
% clear EMGonsite_lock_time;clear EMGonsite_Valid
clear AllBehaviour

%准备延后时间
fs = EEG.srate; % 获取采样率
latency_shift= (232+80) * (fs / 1000); % 结束前232秒一般为应答时间如无法更新


count=length(EEG.event);

for i = 1:1024

    if ~isnan(EMG_onsite(i)) %如有效则更新
    % 新建事件
    new_event.latency = EEG.event(index_trigger10x(i)).latency + EMG_onsite(i)* (fs ); % 在当前触发器之后的1个单位时间加1660
    new_event.type = '55';
    new_event.duration = 0;
    new_event.urevent=count+i;
    
    % 插入到EEG.event中
    EEG.event(end+1) = new_event;%EEG.event(end+1) = EEG.event(i);

    elseif isnan(EMG_onsite(i)) %|| EMG_valid(i)==0 %如无效则按照结束前232秒插入
    new_event.latency = EEG.event(index_trigger10x(i)).latency + EMG_end(i)* (fs )- latency_shift; % 在当前触发器之后的1个单位时间加1660
    new_event.type = '55';
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

index_trigger55 = find(ismember({EEG.event.type}, '55'));
index_trigger66 = find(ismember({EEG.event.type}, '66'));

if length(index_trigger55)~= 1024 && length(index_trigger66)~= 1024%检查数量
    error('length problem, please check');

end
%重新排序
EEG = pop_editeventvals(EEG, 'sort', {'latency', [0]});
%%保存和检查
EEG = eeg_checkset( EEG );
if exp.filter.lowerbound == 0.1
    filename = ['f01' EEG.filename];
elseif exp.filter.lowerbound == 0.01
    filename = ['f' EEG.filename];
end

EEG.filename = filename;

%EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset

%filename = ['aac' exp.name num2str(sub)];
EEG1 = pop_saveset(EEG, filename, exp.filepath);
% Display success message
disp(['RL_File  ' filename '  saved successfully.']);
%% 这个方法总少获取几个 trigger, weifa
% 获取所有事件
% events = EEG.event;
% count=length(EEG.event);

% %找到证据开始点
% index_trigger10x = find(ismember({EEG.event.type}, {'4','44'}));
% index_trigger10 = find(ismember({EEG.event.type}, '10'));
% 
% if length(index_trigger10x)~= 1024 && length(index_trigger10)~= 1024%检查数量
%     error('length problem, please check');
% 
% end
% 
% 
% % 遍历事件查找101-104
% for i = 1:length(events)
%     % 检查是否为trigger 101, 102, 103, 104
%     if ismember(events(i).type, {'4','44'})
%         % 找到最近的trigger 10
%         found_10 = false;
%         for j = i+1:length(events)
%             if strcmp(events(j).type, '10')
%                 found_10 = true;
%                 % 检查101-104和10之间是否有trigger 12或13
%                 found_12_13 = false;
%                 for k = i+1:j-1
%                     if ismember(events(k).type, {'12', '13'})
%                         found_12_13 = true;
%                         % 在trigger 12/13之前1毫秒加trigger 99
%                         % new_event = events(k);
%                         % new_event.type = '55';
%                         % new_event.latency = events(k).latency - 1 * (EEG.srate / 1000); % 前1ms
%                         % events = [events, new_event]; % 插入新的trigger 99
%                         new_event.latency = events(k).latency - 2 * (EEG.srate / 1000); % 
%                         new_event.type = '55';
%                         new_event.duration = 0;
%                         new_event.urevent=count+i;
%                         % 插入到EEG.event中
%                         EEG.event(end+1) = new_event;%EEG.event(end+1) = EEG.event(i);
%                         break;
%                     end
%                 end
% 
%                 % 如果101-104和10之间没有12或13，则在trigger 10之前加trigger 99
%                 if ~found_12_13
%                     % new_event = events(j);
%                     % new_event.type = '55';
%                     % new_event.latency = events(j).latency - 1 * (EEG.srate / 1000); % 前1ms
%                     % events = [events, new_event]; % 插入新的trigger 99
%                     new_event.latency = events(j).latency - 1 * (EEG.srate / 1000); % 前1ms
%                     new_event.type = '55';
%                     new_event.duration = 0;
%                     new_event.urevent=count+i;
%                 end
%                 break;
%             end
%         end
% 
%         % 如果没找到trigger 10，则继续查找下一个101-104
%         % if ~found_10
%         %     continue;
%         % end
%     end
% end
% % 按时间顺序排序事件
% % [~, sorted_idx] = sort([events.latency]);
% % EEG.event = events(sorted_idx);
% index_trigger55 = find(ismember({EEG.event.type}, '55'));
% index_trigger66 = find(ismember({EEG.event.type}, '66'));
% 
% if length(index_trigger55)~= 1024 && length(index_trigger66)~= 1024%检查数量
%     error('length problem, please check');
% 
% end
% %重新排序
% EEG = pop_editeventvals(EEG, 'sort', {'latency', [0]});
% 
% 
% %%保存和检查
% EEG = eeg_checkset( EEG );
% if exp.filter.lowerbound == 0.1
%     filename = ['f01' EEG.filename];
% elseif exp.filter.lowerbound == 0.01
%     filename = ['f' EEG.filename];
% end
% 
% EEG.filename = filename;
% 
% %EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset
% 
% %filename = ['aac' exp.name num2str(sub)];
% EEG1 = pop_saveset(EEG, filename, exp.filepath);
% % Display success message
% disp(['RL_File  ' filename '  saved successfully.']);
end