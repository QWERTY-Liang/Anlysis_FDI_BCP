function [] = TLBF2_cleanEpochs(sub,exp, EEG)

% [EEG1 Indexes] = pop_eegthresh(EEG, 1, 1:128, exp.lowthresh, exp.upthresh, -0.200, 2, 1, 0);
% Interpolate bad channels on a trialo-by-trial basis. This avoids
% rejecting too many trials! 
 
% [EEG aidx] = pop_eegthresh(EEG, 1, 1:exp.nEEGchans, -150 , 150 ,-0.2 , 2, 1, 0);
% EEG.aidx = aidx;
% EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan,-150 , 150 ,-0.2 , 2.0, 1, 0);EEG = pop_TBT(EEG,EEG.reject.rejthreshE,12,0.3,1);
% EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan,-150 , 150 ,-0.2 , 2.0, 1, 0);
% EEG.reject.rejthreshE(129:136,:) = 0;
% EEG = pop_TBT(EEG,EEG.reject.rejthreshE,12,0.3,1);

% EEG.data = EEG.data-EEG.prePulse_baseline;


EEG.chanlocs(129).labels = 'VEOG1';
EEG.chanlocs(130).labels = 'VEOG2';
EEG.chanlocs(131).labels = 'HEOG1';
EEG.chanlocs(132).labels = 'HEOG2';
% 
% EEG1 = pop_eegmaxmin(EEG, [],[], 75, [], 1, 0);
% 
% % EEG1 = pop_TBT(EEG1, EEG1.reject.rejmaxminE , 10, 0.15, 1);
% 
% my_bads = EEG1.reject.rejmaxminE;



%EEG1 = pop_eegthresh(EEG, 1, 1:EEG.nbchan,-75 , 75 ,-1.5 , 1, 1, 0);%这里同一【-1 1】秒查看阈值
EEG1 = pop_eegthresh(EEG, 1, 1:EEG.nbchan,-75 , 75 ,-1 , 0.5,  1, 0);%这里同一【-1 1】秒查看阈值
my_badtrial=EEG1.reject.rejthreshE;%加的定位坏道信息
[EEG1,comrej, badlist]= pop_TBT(EEG1,EEG1.reject.rejthreshE,30,1.0,1);%最后一个变成1可检查
% Max 100% bad epochs per channel --badchan-- Proportion (e.g., 0.3) of max bad epochs per channel.
%Max 20 bad channel per epoch   --badseg - Number of max bad channels per epoch.

%For P1, for example, see correction on  Epoch 14 for ch 95,101
%figure; plot(EEG.data(50,:,14)); hold on ; plot(EEG1.data(50,:,14))
% if ~isempty(badlist.bTrial_num)
%     EEG.data(:,:,badlist.bTrial_num) = []; %delete rejected trials; 
%     EEG.epoch(badlist.bTrial_num) = [];
%     EEG.event = EEG1.event; EEG.urevent = EEG1.urevent;
%     EEG.trials = size(EEG.data,3);
%     EEG.sumData(badlist.bTrial_num,:) = [];
% 
%    % EEG = EP1_makeTrials(sub,exp, EEG,0); %Add structure with trial info for easy sorting; excludes rejected trials.
% 
% else 
%     disp('Debugging')
% end
    
%Rewrite EEG channels on original datafile; keep externals as they were.
%EEG.data(1:128,:,:) = EEG1.data(1:128,:,:);
EEG1.badList_TBT = badlist; %正确答案
EEG1.reject.rejthreshE=my_badtrial;
% Max 0.3 bad epochs per channel --badchan-- Proportion (e.g., 0.3) of max bad epochs per channel.
%Max 20 bad channel per epoch   --badseg - Number of max bad channels per epoch.
%后续加算法提取删掉的trial
EEG1 = eeg_checkset( EEG1 );
filename = ['a' EEG1.filename];
EEG1 = pop_saveset( EEG1, filename, exp.filepath);

 % save as .mat file   
    %save([exp.behpath 'Rejected_' filename ],'my_badtrial','badlist');
end
