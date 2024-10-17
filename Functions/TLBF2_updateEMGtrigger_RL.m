function [EEG1] = TLBF2_updateEMGtrigger_RL(EEG,exp)

% % Load converted data
%EEG = pop_loadset([exp.filepath 'aac' exp.name num2str(sub) '.set']);

% 获取所有事件
events = EEG.event;
count=length(EEG.event);
% 遍历事件查找101-104
for i = 1:length(events)
    % 检查是否为trigger 101, 102, 103, 104
    if ismember(events(i).type, {'101', '102', '103', '104'})
        % 找到最近的trigger 10
        found_10 = false;
        for j = i+1:length(events)
            if strcmp(events(j).type, '10')
                found_10 = true;
                % 检查101-104和10之间是否有trigger 12或13
                found_12_13 = false;
                for k = i+1:j-1
                    if ismember(events(k).type, {'12', '13'})
                        found_12_13 = true;
                        % 在trigger 12/13之前1毫秒加trigger 99
                        % new_event = events(k);
                        % new_event.type = '55';
                        % new_event.latency = events(k).latency - 1 * (EEG.srate / 1000); % 前1ms
                        % events = [events, new_event]; % 插入新的trigger 99
                        new_event.latency = events(k).latency - 1 * (EEG.srate / 1000); % 
                        new_event.type = '55';
                        new_event.duration = 0;
                        new_event.urevent=count+i;
                        % 插入到EEG.event中
                        EEG.event(end+1) = new_event;%EEG.event(end+1) = EEG.event(i);
                        break;
                    end
                end
                
                % 如果101-104和10之间没有12或13，则在trigger 10之前加trigger 99
                if ~found_12_13
                    % new_event = events(j);
                    % new_event.type = '55';
                    % new_event.latency = events(j).latency - 1 * (EEG.srate / 1000); % 前1ms
                    % events = [events, new_event]; % 插入新的trigger 99
                    new_event.latency = events(j).latency - 1 * (EEG.srate / 1000); % 前1ms
                    new_event.type = '55';
                    new_event.duration = 0;
                    new_event.urevent=count+i;
                end
                break;
            end
        end
        
        % 如果没找到trigger 10，则继续查找下一个101-104
        if ~found_10
            continue;
        end
    end
end
% 按时间顺序排序事件
% [~, sorted_idx] = sort([events.latency]);
% EEG.event = events(sorted_idx);

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
end