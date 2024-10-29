function [] = TLBF2_manualChanCheck(sub,exp,EEG1, EEG2)
% This script opens matfiles with single-trial ERP epochs that have been extracted by ExtractEpoch.m,
% which were not yet re-referenced or had any channel interpolated ('raw') but may have been filtered/detrended, 
% and identifies channels with abnormally high (or low!) variance that should be marked for interpolation in the next step 

% NB This script does NOT do the trial-by-trial artifact rejection - that comes later. This step should only 
% flag channels that are bad due to poor connection/ faulty electrode - not ones that are 
% noisy due to blinks/eye movements (those are not bad electrodes - they reflect a bad subject!)

% This procedure is going to seem quite subjective, but it makes you consider each subject/block individually and think through 
% any unexpected glitches you see, and how to salvage as many datasets as possible. Automatic algorithms almost always 
% tend to allow in too little or too much data. A critical principle, though, is to do this in an unbiased way - 
% Our experiments are typically repeated-measures, so it's not a big problem if one subject is noisier than another
% (which is the case anyway whether you use a blanket policy or not)

% We're going to plot the standard deviation of each channel and spot ones that stick out. It's worthwhile to particularly
% scrutinize electrodes that might be particularly important in the final analysis.
exp.ImportantChans = [4 19 85 54 115]; % electrodes around CPz (where CPP is), Fz (CNV) and left and right motor cortex (see 'cap_128_layout_medium.jpg')
% Also mark front-most channels:
exp.FrontChans = [71 72 80 81 93 94 103 70 73 79 82 92 95 102]; % at the very front edge SD can be high not necessarily because channels are bad but because person is blinking /moving eyes. We don't want to interpolate GOOD electrodes that later will HELP us pick up artifacts that we want to throw out.
exp.rerefchan = [19]; % pick a channel to re-reference the data to and get a second picture of variance across channels - important just because sometimes the reference used to load the EEG might itself have been bad during the recording...

% Load in the epoched data
% EEG = pop_loadset([exp.filepath '/P1_fatmc' (exp.name) '_P1.set']);

% Concatenate all epochs and get Standard deviation (SD) per channel
clear SD SDoz
conc = reshape(EEG1.data(1:exp.nEEGchans,:,:),[exp.nEEGchans,size(EEG1.data,2)*size(EEG1.data,3)]); % concatenate all trials

% Now one thing about these raw (as recorded) EEG signals is that the electrodes near the CMS-DRL reference electrodes will tend to be tiny just due to proximity to reference, and even when very noisy you'd miss them in this comparison of S.D.s for that reason.
% So it's worth also looking at SD when referenced to somewhere else - electrode Oz (23) for example - just to make sure you haven't missed any bad channels
conc2 = conc - repmat(conc(exp.rerefchan,:),[exp.nEEGchans,1]);
for q=1:exp.nEEGchans % include externals here as well, because it might be a good idea to check whether those are noisy too
    SD(q,1) = std(conc(q,:))/100; % measure S.D. of each channel
    SD2(q,1) = std(conc2(q,:))/100;
end

% Are there any channels that stick out in terms of standard deviation?
% To check, plot the SD per channel:
figure
subplot(2,1,1); hold on; plot(SD(1:exp.nEEGchans,:));% ylim([0 200]) % we only plot the channels in the cap because external electrodes are often higher variance (e.g. you might be recording EMG) and annoyingly set the scale so you always have to zoom in, and the purpose here is to identify channels for interpolation which is ALWAYS only the 128 cap channels
%title(['subject ' num2str(s) ' ' allsubj{s}])
% mark the important channels specified above - this is just a visual aid to know which you should particularly consider
for e=1:length(exp.ImportantChans)
    plot([1 1]*exp.ImportantChans(e),[0 max(SD(exp.ImportantChans(e),:))],'b'); 
end
% mark the front edge channels as well, again a visual aid to know where blinks are likely to create higher variance (through no fault of the electrodes)
for e=1:length(exp.FrontChans)
    plot([1 1]*exp.FrontChans(e),[0 max(SD(exp.FrontChans(e),:))],'k'); % front edge channels in BLACK
end

%Identify candidate bad channels
candidates = find(SD(1:exp.nEEGchans,:) > 75 | (SD(1:exp.nEEGchans,:) < 1 & (SD2(1:exp.nEEGchans,:) < 1)));
exclCandidates = setdiff(candidates,exp.FrontChans);

for e=1:length(exclCandidates)
    plot([1 1]*exclCandidates(e),[0 max(SD(exclCandidates(e),:))],'r'); % candidate noisy channels in RED
end

subplot(2,1,2); hold on; plot(SD2(1:exp.nEEGchans,:)); %ylim([0 200])
% mark the important channels specified above - this is just a visual aid to know which you should particularly consider
for e=1:length(exp.ImportantChans)
    plot([1 1]*exp.ImportantChans(e),[0 max(SD(exp.ImportantChans(e),:))],'b'); 
end
for e=1:length(exp.FrontChans)
    plot([1 1]*exp.FrontChans(e),[0 max(SD(exp.FrontChans(e),:))],'k'); % front edge channels in BLACK
end

ch2interp = exclCandidates; % makes an empty cell array for filling in the bad channels for this subject

disp(['Channels to interpolate: ' num2str(exclCandidates')]);
prompt = 'Confirm interpolation selection? y/n [y]: ';
str = input(prompt,'s');


if strcmp(str,'y')
    %Interpolate in SL epochs
    EEG1 = eeg_interp(EEG1,ch2interp);
    EEG1.interpolated = ch2interp;
    EEG = eeg_checkset( EEG1 );
    filename = ['i' EEG1.filename];
    EEG1 = pop_saveset( EEG1, filename, exp.filepath);
    clear EEG1;
     
    %Interpolate same channels in RL epochs
    EEG2 = eeg_interp(EEG2,ch2interp);
    EEG2.interpolated = ch2interp;
    EEG = eeg_checkset( EEG2 );
    filename = ['i' EEG2.filename];
    EEG2 = pop_saveset( EEG2, filename, exp.filepath);

else
    disp('Interpolation aborted');
end


%% Another way to see if a channel that seems bad for certain blocks are actually bad channels, or is it a 
%% bad subject doing something whacky on some trials, is to get SD per channel per trial:
SDct = squeeze(std(EEG.data(1:exp.nEEGchans,:,:),[],2)); % after squeezing, it's the SD of each trial for 128 channels x numtrials
figure; imagesc(SDct,[0 35]);   % the last input of imagesc sets the color axis limits, and I chose them based on typical SD values from the first plot,
% % Vertical stripes of high SD mean lots of channels are bad for a limited set of trials (so not a case for interpolation - just let artifact rejection criterion later exclude those trials) 
% % Horizontal stripes of high SD means lots of trials having high SD for those specific channels - that's a
% % case for interpolation

return; % this is here because for a given subject you want to run the below lines individually
