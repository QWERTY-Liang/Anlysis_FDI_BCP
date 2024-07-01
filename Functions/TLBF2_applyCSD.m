function [] = TLBF2_applyCSD(sub,exp, baseline, epoch)

%Relies on CSD toolbox
%See https://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/tutorial.html

% if strcmp(epoch, 'P1')
% EEG = pop_loadset([exp.filepath '\ab_cICAtrial_eyeB_rej_riP1_f01atmc' exp.name '_P' num2str(sub) '.set']);
% elseif strcmp(epoch, 'RL')
%     EEG = pop_loadset([exp.filepath '\b_cICAtrial_b_RL_f01riatmc' exp.name '_P' num2str(sub) '.set']);
% 
% end

%Load CSD transform
load('CSD_coords_biosemi128.mat');
% TO DO: check that these coordinates are obtained using the right lambda
% and m parameters - see tutorial. 

% %% Baseline prePulse and plot all channels
% figure; plot(mean(EEG.data(1:128,:,:),3)');
% 
% if ~contains(EEG.filename, 'cICA')
%     if strcmp(baseline, 'preStim')
%         for e = 1:length(EEG.epoch)
%             EEG.data(:,:,e) = EEG.data(:,:,e) - EEG.preStim_baseline(:,e);
%         end
%     elseif strcmp(baseline, 'prePulse')
%         for e = 1:length(EEG.epoch)
%             EEG.data(:,:,e) = EEG.data(:,:,e) - EEG.prePulse_baseline(:,e);
%         end
%     elseif strcmp(baseline, 'action')
%         for e = 1:length(EEG.epoch)
%             EEG.data(:,:,e) = EEG.data(:,:,e) - EEG.actionBaseline(:,e);
%         end
%     end
% else
%     if strcmp(baseline, 'preStim')
%         for e = 1:length(EEG.epoch)
%             EEG.data(:,:,e) = EEG.data(:,:,e) - EEG.preStim_baseline_ICA(:,e);
%         end
%     elseif strcmp(baseline, 'prePulse')
%         for e = 1:length(EEG.epoch)
%             EEG.data(:,:,e) = EEG.data(:,:,e) - EEG.prePulse_baseline_ICA(:,e);
%         end
%     elseif strcmp(baseline, 'action')
%         for e = 1:length(EEG.epoch)
%             EEG.data(:,:,e) = EEG.data(:,:,e) - EEG.actionBaseline_ICA(:,e);
%         end
%     end
% end
        
figure; plot(mean(EEG.data(1:128,:,:),3)');
EEG = eeg_checkset( EEG );

%% Apply CSD to epoched, BASELINED data, but only to EEG channels
X = [];
tic();
for e = 1:length(EEG.epoch)
    X(:,:,e) = CSD(EEG.data(1:128,:,e), G,H);
end
toc();
%This takes about 3 minutes per participant!
EEG.data(1:128,:,:) = X;

EEG = eeg_checkset( EEG );
filename = ['csd_' baseline EEG.filename];
EEG = pop_saveset( EEG, filename, exp.filepath);



end