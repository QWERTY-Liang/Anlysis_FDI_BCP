%%  6.1 EMG analysis
% Author: Liang Tong
% Date: 12/7/2024 

%to update list
%1.function all the code to make main script simplier
%% Toolbox requirements: 
clc
clear all
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis'));% current folder

%% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);


%% 1. Load rl_BB pre-cue-188
% col 1: subject num
% col 2: block num
% col 3: contrast (high or low contrast)
% col 4: muscle (1=FDI, 2=BCP)
% col 5: trial outcome (1=correct, 2=error, 4=no response, 3 too early, 6 slow, 5 wrong muscle)
% col 6: RT in sec
% col 7: evshowtime =[]   % 0.2sec additional time after response
% col 8: participant response, 1 = left, 2 = right 0= on response
% col 9: correct response, 1 = left, 2 = right
% col 10: first tilt used for SSVEP
%证据前0.188秒 
% cue前0.188秒 col:12
% 整段平均 col:13
% response 前0.188秒 col:14
addpath(genpath(exp.finalpath));

load TL_ALL_include_EMG_updated.mat% EMG missing BCP totEMG

load TL_ALL_include_EMG% in EMGfolder

load([exp.behpath exp.name 'EMG_brusts'])
%%

EMGonsite=9999*ones(length(AllBehaviour_new),1);
changemark=zeros(length(AllBehaviour_new),1);

for trial=1:6144
    rt=AllBehaviour_new(trial,6);% matlab 记录的反应时间
    %FDIorBCP=AllBehaviour_new(trial,4); %1为用了FDI, 2 为BCP

    if AllBehaviour_new(trial,8)==2%右手回答的情况
        % Extract the start times for the right hand in the current trial
        start_times_right = movement_start_times_right{trial};

        % Find all start times that are less than rt
        valid_times = start_times_right(start_times_right < rt);
        if movement_counts_right(trial)>0

            % Check if there are any valid times
            if ~isempty(valid_times)
                % Find the time closest to rt
                [~, idx] = min(abs(valid_times - rt));
                if (rt-valid_times(idx))<0.4
                EMGonsite(trial) = valid_times(idx); % Store the closest valid time
                changemark(trial)=1;%记录更新
                else
                % If no valid time is found, set EMGonsite to rt
                EMGonsite(trial) = rt;
                changemark(trial)=-1;%记录更新
                end
            else
                % If no valid time is found, set EMGonsite to rt
                EMGonsite(trial) = rt;
                changemark(trial)=-1;%记录更新
            end

        else
EMGonsite(trial) = rt;
        end

    elseif AllBehaviour_new(trial,8)==1%左手回答的情况，信号质量差

        % Extract the start times for the right hand in the current trial
        start_times_left = movement_start_times_left{trial};

        % Find all start times that are less than rt
        valid_times = start_times_left(start_times_left < rt);
        if movement_counts_left(trial)>0

            % Check if there are any valid times
            if ~isempty(valid_times)
                % Find the time closest to rt
                [~, idx] = min(abs(valid_times - rt));
                if (rt-valid_times(idx))<0.4
                EMGonsite(trial) = valid_times(idx); % Store the closest valid time
                changemark(trial)=1;%记录更新
                else
                % If no valid time is found, set EMGonsite to rt
                EMGonsite(trial) = rt;
                changemark(trial)=-1;%记录更新 
                end
            else
                % If no valid time is found, set EMGonsite to rt
                EMGonsite(trial) = rt;
                changemark(trial)=-1;%记录更新
            end
        else
            EMGonsite(trial) = rt;

        end


    end

end



%new EMG onsite data created
%col 11 是更新后的 EMG onsite time
%col 12 是与原来时间的提前量
%col 13 是左手brust计数
%col 14 是右手brust计数
%col 15 是是否更新

Tdiff=AllBehaviour_new(:,6)-EMGonsite;
AllBehaviour_EMGonsite=[AllBehaviour,EMGonsite,Tdiff,movement_counts_left,movement_counts_right,changemark];

save([exp.behpath exp.name 'EMG_brusts_update'], 'AllBehaviour_EMGonsite');