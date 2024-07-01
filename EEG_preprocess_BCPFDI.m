%%--------------------------------------------------------------------------------------------
%to-update list
%1. in annotate data, update trial ID and Meta data

%%  2. EEG analysis script
% Author: Liang Tong
% Date: 4/6/2024 

%% Toolbox requirements: 
%Biosig: for reading BDF file

clc
clear
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis'));% current folder

%% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);

%% Run analysis
eeglab

%% STEP 0: Convert data to EEGlab format
% Convert data from bdf to .set and .fdt
for sub = exp.sub_id(1:end)

    TLBF2_convertData(sub, exp);

end

%% STEP 1: Append all data
for sub = exp.sub_id(1:end)

    TLBF2_appendData(sub, exp);

end

%% STEP 2: Annotate dtata (update later)
%to update later
% 1. Add behavioural data matrices with information about each trial and event
% to the EEG structure.
for sub = exp.sub_id(1:end)

    TLBF2_annotateData(sub, exp);

end

%% STEP 3: Filter data
% Filter data using the 
for sub = exp.sub_id(1:end)
    EEG = TLBF2_filter(sub,exp);
end

%% STEP 4: Epoch, detect bad channels & interpolate
% The manualChanCheck function will plot the SD for all channels with the
% online reference (upper plot) and re-referenced to POz (lower plot) just
% in case the original reference was noisy. 

% The vertical black lines indicate channels that are likely to show high SDs due to
% eyeblinks (these include the two most frontal rows of electrodes). 

% The blue lines indicate particularly important channels that we may need
% for analysis. These are crucial, and its important to note if they are
% noisy. 

% The red lines are the channels that will be suggested for interpolation
% to you. These will be noisy channels (SD>50 across the whole experiment) 
% and that are NOT frontal ones. Take a look at the data and decide whether
% you agree. If yes, type 'y' on the console to accept interpolation. 
% If no, type 'n'.

% Please keep a log of which channels you reject for each
% participant.

for sub = exp.sub_id(1:end)
    [EEG1, EEG2] = TLBF2_epochData(sub,exp);  %Epoch around evidence and around response
    TLBF2_manualChanCheck(sub,exp,EEG1,EEG2);
    clear EEG1; clear EEG2; close all;
end
%logbook for bad channel (100uV thresh)
%TL1: 22  23  27  28  37  40
%TL2: /
%TL3: /
%TL4: /
%TL5: /
%TL6: /


%% STEP 5: Reject epochs with eyeblinks during pulses, reject artefacts, & compute baselines & create trial list
%to update list: add metadata and trial sorting
% now only re-refernce
% In this step, and after we have interpolated bad channels, we
% re-reference the data to the common average (i.e. average of all
% electrodes). 

% Then, we reject trials where participants blinked during one of the two
% pulses. Trials are rejected if two criteria are met: 
% 1) The VEOG shows fluctuations > 150 uV during the pulses, AND 2) artefacts
% of > 80 uV are found in any of the frontal channels. Criterion 2) is
% added because VEOG channels are sometimes noisy, and that criterion would
% lead to the rejection of too many trials. 

% Finally, we create a trial structure that contains a list of all relevant
% experimental factors we may want to group trials by. 

for e =1:2
    epoch = exp.epochs{e}
    allRej = [];
    for sub = exp.sub_id(1:end)
        EEG = TLBF2_reReference(sub,exp,epoch);
        % if strcmp(epoch,'SL')
        %     [EEG allRej] = EP1_removeEyeblinkEpochs(sub,exp,EEG,allRej); %Reject epochs with eyeblinks during pulses
        % end
        % EEG = TLBF2_makeTrials(sub,exp, EEG,1); %Add structure with trial info for easy sorting; excludes rejected trials.
    EEG = pop_saveset(EEG, [EEG.filename], exp.filepath);
    end
    %save('GapsTask_eyeBlink_rej', 'allRej');
end


%% STEP 6: run ICA and remove eye movement components 
% In this step, we run an Independent Component Analysis to identify eye
% movements (blinks and saccades).

% We run the ICA on downsampled data to speed up the process, and we then 
% save the ICA components to the original data with the full sampling rate.
% To visualise the components, open the eeglab UI (type eeglab on the
% matlab console), and load the file called ICA_xxx for your participant. 
% Then click on Plot --> Component maps --> In 2-D, and wait for the 130
% components to be plotted. 

% Then, see if you can identify which components correspond to eyeblinks
% and which ones to saccades. 
% Vertical eye movements typically have a frontal topography, whereas
% horizontal eye movements (saccades) show a left-to-right topography.
% These are easy to identify once you get used to it! 
% Here is a useful guide to understanding ICA: 
% https://eeglab.org/tutorials/06_RejectArtifacts/RunICA.html 

% Once you have two candidate components for the blinks and saccades, look
% at what happens if you remove them. 
% On the EEGLAB UI, click Tools --> Remove components from data. 
% On the panel, type the numbers of components that you want to reject. 
% On the pop-up panel, click "Plot single trials" to visualise what happens
% to the data after removing your components. If you picked the right
% components, you should see that the red line is flat where there are
% eyeblink artefacts and horizontal eye movements. The red line should
% overlap almost perfectly with the blue one at all other times. It is VERY
% important to always check this step, becaues sometimes the ICA algorithm
% will not work perfectly and may introduce noise when removing components.
% 

% Keep a log of the components you identify for each participant (possibly
% take screenshots as well) and make sur eyou note down how clean the
% artefact removal is. 

% Once you are done and you think you have identified the right components,
% click "Accept" on the pop up menu and save your corrected file by adding
% a 'c' at the beginning of your original file (e.g. cICAtrial_eyeB...). 

for e = 1:2
    epoch = exp.epochs{e}
    for sub = exp.sub_id(1:end)
        % if strcmp(epoch,'SL')
        %     EEG = pop_loadset([exp.filepath 'ri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
        % elseif strcmp(epoch, 'RL')
            EEG = pop_loadset([exp.filepath 'ri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
        % end
        TLBF2_runICA(sub,exp, EEG);   
    end
end

%% STEP 6.5: (optional):  remove eye movement components (using ICA lable plugin)

% bug fixed: last 4 channel dropped (ICAlable evalc 被注释掉)


for e = 1:2
    epoch = exp.epochs{e}
    for sub = exp.sub_id(1:end)
        % if strcmp(epoch,'SL')
        %     EEG = pop_loadset([exp.filepath 'ri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
        % elseif strcmp(epoch, 'RL')
            EEG = pop_loadset([exp.filepath 'ICAri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
        % end
        TLBF2_lableICA(sub,exp,EEG);
 
    end
end
%% STEP 7: reject all epochs where there are fluctuations >150uV in any channel
%这一步删掉了一些trial，需后期定位把trial的idx找出来（可从TBT里找例子）
%这里统一用了0点附近-1.5-1秒区域
% This is almost the last step! We will now baseline the data before the
% first pulse, and then do a final cleaning step. 

% For this cleaning step, we are using a toolbox that allows us to
% interpolate bad channels on a trial-by-trial basis. 
% The toolbox can cause some problems occasionally

% At the moment, it will ask you to accept the interpolation & trial
% rejection step. It is a good idea for you to plot the EEG data and the
% matrices for each participant to see what is being rejected. 

% It will normally always reject the 8 external channels, because they have
% not been baselined and are always out of bounds - this is OK. We are 
% keeping the original data. If other channels are being rejected as well, 
% keep a note of them. 


for e = 1:2
    epoch = exp.epochs{e};
    for sub = exp.sub_id(1:end)
        % if strcmp(epoch,'SL')
        %     EEG = pop_loadset([exp.filepath 'ri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
        % elseif strcmp(epoch, 'RL')
            EEG = pop_loadset([exp.filepath 'cICAri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
        % end
         [EEG1 EEG2] = TLBF2_baselineData(sub,exp, EEG, epoch); %Calculate baselines and save in a single data file; by default, subtract pre-Pulse/stlmu for 'SL' per-action for RL
       %tbt2里面改了，129-136全为坏电极默认，这步后无外接电极
       %每个人1024有8个trial左右扔掉，可以更严格
         TLBF2_cleanEpochs(sub,exp,EEG1); % Threshold; remove epochs with > 150uV drifts in any channel
       TLBF2_cleanEpochs(sub,exp,EEG2);
    end
end

%% STEP 8: CSD filtering 
% Relying on CSD toolbox (https://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/)
% Final step: Apply CSD transformation to the data 
for e = 1:2
    epoch = exp.epochs{e};
    for sub = exp.sub_id(1:end)

         EEG = pop_loadset([exp.filepath 'ab_cICAri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
        TLBF2_applyCSD(sub,exp, EEG)
         EEG1 = pop_loadset([exp.filepath 'abb_cICAri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
        TLBF2_applyCSD(sub,exp, EEG1)

    end
end

%% STEP 9: Morlet - wavelet transform 


