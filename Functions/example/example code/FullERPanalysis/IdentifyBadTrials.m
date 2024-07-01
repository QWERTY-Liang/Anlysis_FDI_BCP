% This script now goes through each subject and figures out the best way to identify 'bad' trials to be excluded from averaging

clear all
allsubj = {'P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16'}; % remember there are 8 more!

bigmatfolder = 'bigmats/';

% load up and define some of the same things we need:
load chanlocsBioSemi128;
nchan = 128; % number of channels
next = 8;  % number of externals
fs=1024; % sample rate
VEOGchans = [129 130]; %  CHECK THAT THESE ARE INDEED YOUR VEOG CHANNELS [upr lwr]

% Before sorting trials and averaging by condition, we need to figure out what analysis time window and artifact
% rejection (AR) criteria will give us the best data quality across as many subjects as possible.
% AR criteria by our very simple threshold method consist simply of a time window in which to check for
% artifacts, and a threshold level that marks a trial as having an artifact if exceeded within that window.
% The first thing to consider is blink behaviour because it's the biggest artifact source. If you instructed
% the subjects well, they will have waited until the stimulus has fully ended or until a good 0.5-1 sec after
% responding before they get their blink in, but sometimes they blink very quickly and if you set too wide an
% AR window in this case, loads of trials will be rejected even if the time between stimulus and response is
% mostly clean. When there are several such subjects in a dataset, it can be best to artifact-reject just
% based on a relatively narrow window from stimulus onset and a narrow window just prior to response - this is what I arrive
% at below in fact. It is far better to instruct your subjects well than to have to scrape around for trials by
% adjusting windows like this! 
% A standard approach is to reject a trial if there is ANY sample within the WHOLE epoch window on ANY electrode that exceeds a threshold
% value. Other approaches are to consider only a subset of electrodes that are particularly relevant, or as above, check a narrower time
% window WITHIN the overall epoch. A compromise I've used before is to set a lower threshold on the most relevant electrodes to make sure
% they are artifact free, and a higher threshold on all others. There is no one absolutely correct way, and sometimes adjusting these things
% makes little difference. But sometimes it makes a substantial difference so it's important to test and explore your artifact rejection
% settings.

% If whatever window/electrodes you use results in too much rejection (typically >30% considered too much, but it really depends on number of trials overall), 
% you should first check if there are particular electrodes that are the culprit most of the time. If they are isolated and stick out among
% neighbours, maybe that's actually a bad electrode that you should have interpolated, so you can go back and add to the ch2interp list and
% re-interpolate.

% %%
% % We will now 'draft' some artifact detection settings and check:
% % 1) how many trials get rejected? 
% % 2) Is it because of blinks or other artifacts, and if blinks, can we adjust the epoch window of our analysis to avoid the blinks?
% % 3) If a lot of rejection is due to non-blink artifacts, are there particular channels mostly responsible and should we
% % go back to previous step and add them in ch2interp for that subject and re-run their interpolation?
% 
% % Here are reasonable artifact-rejection settings to begin with - thresholds, channels to check and time window to check:
% blinkTh = 100; % threshold for detecting a blink in VEOG upper-lower difference signal. 100 uV usually works well
% artifTh = 80; % this is for SCALP channels - note that we are doing this before CSD transformation. 80 is typical but might need adjusting
% ARchecklim = [-100 2000]; % in msec, the limits of the time window relative to the event in which you will check for an artifact
% % Here I'm trying a relatively wide window that should include most RTs, and we'll see if we need to tighten it
% ARchans = [1:128]; % Artifact rejection channels - only these will be checked for exceeding threshold. All 1:128, or just some key ones?
% % Usually it's best to use all channels because we have to plot topographies and don't want artifacts on certain
% % channels muddying the topographies, but if the data are poor quality and hard to salvage, we might only check the
% % channels where key signals are expected to show up, e.g. CPP, beta, and could set higher thresholds on the rest.
% 
% % I first run the below loop with the above settings and assess the situation, then I try alternative settings in the
% % next cell.
% 
% blinkTms_all = []; blinkTms_allR = []; % (initialisation) these will contain the times relative to stimulus and to response ('R') at which blinks occurred on the trials
% cumMaxAbsFig = figure; hold on % we'll put the cumulative distribution of the maximum absolute values per trial in one figure for all subjects together
% for s=1:length(allsubj)
% 
%     subjID = allsubj{s};
%     disp([subjID '...'])
% 
%     load([bigmatfolder subjID '_intp']) % the interpolated, but not yet CSD'd data
% 
%     % We'll put some things in a big figure for each subject to summarise their artifact tendencies:
%     figure; set(gcf,'Position',[80 80 1000 700]);
% 
%     % Now the data are still referenced to the original online reference on the head. That means channels close to that
%     % will be potentially small by proximity only, so we will average-reference the data to evaluate artifacts on each
%     % channel with respect to the overall mean across electrodes:
%     erp(1:nchan,:,:) = erp(1:nchan,:,:) - repmat(mean(erp(1:nchan,:,:)),[nchan,1,1]);
% 
%     % First analyse blinks. Let's assume most of the time window has eyes-open so the median value will
%     % represent that state, to be used as a baseline:
%     VEOG = squeeze(erp(VEOGchans(1),:,:)-erp(VEOGchans(2),:,:)); % Upper elec minus lower elec; so now it's time x trials
%     VEOG = VEOG - repmat(median(VEOG),[size(VEOG,1),1]); % "subtract the value estimated as the eyes open state, so blinks are seen as deviations from this
%     % now plot VEOG to investigate blinking behaviour:
%     subplot(2,3,1); plot(t,VEOG); hold on; plot(t,ones(size(t))*blinkTh,'k'); title(['Subject ' allsubj{s}]) % plot the VEOG waveforms for all trials. e.g. you can see P01 makes lots of blinks but quite seldom during the key time period just following a target onset.
%     % For each trial get the time at which the blink threshold is exceeded, to look at the distribution of
%     % times subjects tend to blink:
%     blinkTms=[]; 
%     for n=1:size(VEOG,2) % for each single trial
%         blnk=find(abs(VEOG(:,n))>blinkTh,1); % this captures the FIRST blink occurring in the epoch window - keep in mind there might be more after that. This stll provides a good indication
%         if isempty(blnk), blinkTms(n)=nan; else, blinkTms(n)=t(blnk); end; 
%     end; 
%     subplot(2,3,2); hist(blinkTms,[t(1):100:t(end)]); title(['S-lock']); 
%     subplot(2,3,3); trl = find(RT>0 & RT<2); hist(blinkTms(trl)-RT(trl)*1000,[-1500:100:1500]); title('R-lock') % useful to see also when the blinks are happening relative to RT. Just including trials with a reasonable RT of <2 sec here
%     blinkTms_all = [blinkTms_all blinkTms]; % and append to a big vector of blink times to check across all subjects
%     blinkTms_allR = [blinkTms_allR blinkTms(trl)-RT(trl)*1000];
%     % Look at subject P01 for example: they seem to blink a lot but only just after they pressed the button!
% 
%     % the following will count up all of the trials with artifacts detected, whether blinks or other large activity on
%     % cap channels. First get the window:
%     artcheckwin = find(t>=ARchecklim(1) & t<=ARchecklim(2)); % indices of the timepoints for the artifact check window whose limits (min and max) are defined above
%     % now check for blinks in that window:
%     blink = max(abs(VEOG(artcheckwin,:))) > blinkTh; % this will give a logical 1 if this threshold is exceeded on a certain trial, 0 otherwise
%     % and then artifacts on any of the channels specified in ARchans above:
%     maxabs = squeeze(max(abs(erp(ARchans,artcheckwin,:)),[],2)); % This computes max absolute value per electrode per trial (dimensions: elec x trial)
%     art = maxabs-artifTh > 0;  % Across trials flag art = 1 if any timepoint exceeds artif thresh, 0 if fine (elec x trial)
%     artifact = sum(art,1); % How many electrodes had an artifact on each TRIAL (1d vector, trials). So this will be 0 if no artifact on any and can be used as a logic true/false
%     numartPerChan = sum(art,2); % How many trials had an artifact for each ELECTRODE
%     subplot(2,3,4); bar(numartPerChan); title(['number of artifact trials per channel']); grid on
%     numartPerChan1 = sum(art(:,find(~blink)),2); % How many trials THAT HAD NO BLINK had an artifact for each ELECTRODE
%     subplot(2,3,5); bar(numartPerChan1); title(['On the trials with no blinks:']); grid on
%     % the last plot is important for seeing if there are channels that we need to go back and interpolate in the
%     % previous step. If channels are exceeding artifact thresholds just because of blinks, then they're not bad channels.
%     % In P01 for example, the artifact trials per channel (bottom left) looks spiky, showing certain electrodes exceed threshold far more than the rest, but
%     % these are actually the ones right at the eyes, e.g. C16 (80) and C29 (93), and when you only plot trials with no blink detected at
%     % VEOG (bomttom middle), it shows there is no obvious noisy electrode that isn't just capturing blinks.
% 
%     % For each subject, proportion of trials that have blinks: 
%     propBlink(s) = length(find(blink))/length(RT);
%     % proportion of trials that have artifacts on other channels (bearing in mind this included frontal channels also picking up blinks)
%     propArtifact(s) = length(find(artifact))/length(RT);
%     % proportion of trials that would be kept if these were the artifact/blink detection criteria:
%     propKept(s) = length(find(~blink & ~artifact))/length(RT); 
% 
%     % Now looking at just non-blink trials, how appropriate is the artifact threshold on the scalp channels, in terms of striking the
%     % balance between retaining a good number, but kicking out the worst trials?
%     figure(cumMaxAbsFig);
%     nonblinktr = find(~blink); n_nbt = length(nonblinktr);
%     maxmaxabs = max(maxabs(:,nonblinktr)); % get the maximum across electrodes of all the maximum absolute values for each trial
%     plot([1:n_nbt]/n_nbt,sort(maxmaxabs)); % sort from best trial to worst trial and plot (only the non blink trials - this is to assess choice of artifTh)
% end
% % spruce up this figure, which orders the maxmaxabs across trials and gives a sense of the proportion of trials that
% % would make the cut, AND how that might potentially change if you were to adjust the threshold
% figure(cumMaxAbsFig); title('ordered maxmaxabs for each subject'); legend(allsubj(1:s)); xl=xlim; plot(xl,[1 1]*artifTh,'k--') % will need to zoom in!
% ylabel('maxmaxabs'); xlabel('proportion of trials')
% figure; set(gcf,'Position',[100 180 1000 400]); subplot(1,3,1), bar(propBlink); title('propBlink'),  subplot(1,3,2), bar(propArtifact); title('propArtifact'); subplot(1,3,3), bar(propKept); title('propKept'); 

% %% assessment:
% % First of all, looking at the bar plots of proportions per subject, "proportion kept" (not rejected based on either a blink or a large
% % amplitude on any cap electrode) using this long artifact check window (ARchecklim) ranges from excellent (97%) to terrible (15%). 
% % You'd really prefer most proportions kept to be up higher than 80%.
% % As you can see from the individual-subject VEOG plots, some subjects are good, some are terrible blinkers. For many subjects you can see that they
% % rarely blink during the first 0.5-1 sec following evidence onset but relative to the response (rightmost plot), they blink very soon afterwards. 
% % Let's check the histograms of blink times across all subjects (pooled): 
% figure; subplot(1,2,1), hist(blinkTms_all,[t(1):100:t(end)]); title('S-lock'); subplot(1,2,2); hist(blinkTms_allR,[-2500:100:1500]); xlim([-2500 1500]); title('R-lock')
% 
% % Subjects look like they are trying to time a blink after their response, but often, too soon and sometimes ON the response. 
% % The experimenter should have instructed the subjects more clearly to say "try to wait at least half a
% % second after your button click before you blink, if you need to." 
% % With this dataset, we will clearly have to carefully tailor the settings so as to squeeze as many trials as possible out of each subject.
% % The goal is always to retain as many trials as possible but to avoid letting in noisy trials, to get the most reliable average
% % The first thing we should try is restricting the artifact check window, since obviously a lot of blinks are happening within 2000 ms of
% % the target onset.
% 
% % more notes on individuals:
% % You sometimes have to zoom into the VEOG because some craziness on a minority of trials is making it zoom out a lot , e.g. P15. By and
% % large, these figures allow us to confirm that the VEOG signals were mostly recorded well and doing their job.
% % P11 has their VEOG electrodes upside down. It's fine since we check absolute value - this is why we do that!
% % In most cases there is not all that much additional artifact rejection beyond the blinks, which is good.


%%

% This is our plan:
% We will choose target-locked and response-locked windows that are just wide enough to capture the dynamics 
% of the decision process just after target onset and just before response, and only artifact-reject based on data in those narrow windows

epoch_limits_msTG = [-53 747];    % Target-locked epoch; again using integer number of SSVEP cycles (1000/18.75 = 53.33)
tts = round(epoch_limits_msTG(1)/1000*fs):round(epoch_limits_msTG(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
tt = tts*1000/fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds
epoch_limits_msR = [-427 53];    % Response-locked epoch; again using integer number of SSVEP cycles (1000/18.75 = 53.33)
trs = round(epoch_limits_msR(1)/1000*fs):round(epoch_limits_msR(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
tr = trs*1000/fs; % hence timebase in milliseconds
% Note that in Response-locked waveforms time '0' is the time of response, and we usually don't care so much about 
% activity after that time. Our strategy is to cut out the shorter T-locked and R-locked epochs from EACH of the 
% single stimulus-locked epochs. We have to do it on a trial by trial basis (inside the loop below) because RT is different from trial to trial!

% And the thresholds and channels to check (let's keep them the same as above for now):
blinkTh = 100; % VEOG blink threshold 
artifTh = 80; % artifact threshold for SCALP channels 
ARchecklim = [-100 2000]; % in msec, defines artifact check window
ARchans = [1:128]; % Artifact rejection channels

cumMaxAbsFig1 = figure; hold on % we'll put the cumulative distribution of the maximum absolute values per trial in one figure for all subjects together
clear avERP avERPr
for s=1:length(allsubj)
    
    subjID = allsubj{s};
    disp([subjID '...'])
    
    load([bigmatfolder subjID '_intp'])

    figure; set(gcf,'Position',[80 80 700 400]);

    % Extract response-locked ERPs, create a matrix of single-trial erps similar to 'erp' but with each single trial time-locked to the response 
    ntr = length(RT); % number of trials, which we can get from the length of the RT vector
    erpr = nan(size(erp,1),length(tr),ntr); % initialise to nan - this is the correct matrix size and we'll fill it in in the single-trial lop below
    validrlock = zeros(1,ntr); % make a vector that says whether each trial is a valid trial for response locking (1) or not (0)
    for n=1:ntr
        [blah,RTsamp] = min(abs(t-RT(n)*1000)); % find the sample point closest to the RT for this trial. Remember RT is in sec
        if RTsamp+trs(1) >0 & RTsamp+trs(end)<=length(t) & ~isnan(RT(n)) % a trial is only validrlock=1 if the RT lies within the erp epoch you're using (sometimes there's no response at all, and the RT for these misses might show up as 4.5 sce or something, which is actually the response to the following trial! or you might decide to code missed RTs as NaN. This catches all those
            erpr(:,:,n) = erp(:,RTsamp+trs,n); % Extract response-locked epoch
            validrlock(n)=1;
        end
    end
    % note we don't have to do the same for target locked because the long 'erp' epoch is already target
    % locked - all we have to do is only artifact-check in the narrower range, and plot that range.
    
    % find artifact trials:
    artcheckwinTG = find(t>=epoch_limits_msTG(1) & t<=epoch_limits_msTG(2)); % indices of the timepoints for the artifact check window
    % First blinks:
    VEOG = squeeze(erp(VEOGchans(1),artcheckwinTG,:)-erp(VEOGchans(2),artcheckwinTG,:)); % Upper elec minus lower elec; so now it's time x trials
    % not going to bother with baseline-correcting VEOG here because we will check absolute deviation,+ or - value greater than threshold
    VEOGr = squeeze(erpr(VEOGchans(1),:,:)-erpr(VEOGchans(2),:,:)); % response-locked VEOG, full timeframe
    blink = max(abs([VEOG;VEOGr])) > blinkTh; % this will give a logical 1 if this threshold is exceeded for EITHER the T-locked or R-locked VEOG, 0 otherwise
    % and then artifacts on any of the channels specified in ARchans above
    maxabs = squeeze(max(abs([erp(ARchans,artcheckwinTG,:) erpr(ARchans,:,:)]),[],2)); % This computes max absolute value per electrode per trial (dimensions: elec x trial)
    art = squeeze(max(abs([erp(ARchans,artcheckwinTG,:) erpr(ARchans,:,:)]),[],2)-artifTh > 0);  % art is elec x trial, 1 if any timepoint exceeds artif thresh, 0 if fine. Again, concatenated t-lock and R-lock so will detect artifact in EITHER
    numart = sum(art,2); % How many trials had an artifact for each ELECTRODE    
    artifact = sum(art,1); % How many electrodes had an artifact on each TRIAL
    numartPerChan = sum(art,2); % How many trials had an artifact for each ELECTRODE
    subplot(1,2,1); bar(numartPerChan); title([subjID ': # artifact trials per channel']); grid on
    numartPerChan1 = sum(art(:,find(~blink)),2); % How many trials THAT HAD NO BLINK had an artifact for each ELECTRODE
    subplot(1,2,2); bar(numartPerChan1); title(['On the trials with no blinks:']); grid on
      
    % For each subject, proportion of trials that have blinks: 
    propBlink(s) = length(find(blink))/length(RT);
    % proportion of trials that have artifacts on other channels (bearing in mind this included frontal channels also picking up blinks)
    propArtifact(s) = length(find(artifact))/length(RT);
    % proportion of trials that would be kept if these were the artifact/blink detection criteria:
    propKept(s) = length(find(~blink & ~artifact))/length(RT); 
        
    % Now looking at just non-blink trials, how appropriate is the artifact threshold on the scalp channels, in terms of striking the
    % balance between retaining a good number, but kicking out the worst trials?
    figure(cumMaxAbsFig1);
    nonblinktr = find(~blink); n_nbt = length(nonblinktr);
    maxmaxabs = max(maxabs(:,nonblinktr)); % get the maximum across electrodes of all the maximum absolute values for each trial
    plot([1:n_nbt]/n_nbt,sort(maxmaxabs)); % sort from best trial to worst trial and plot (only the non blink trials - this is to assess choice of artifTh)

    % If we're happy with the results of all the above, re-run this with the next bit uncommented, to save the artifact/blink info:
    % add details to anapar, the analysis parameters record:
    anapar.tt = tt; anapar.epoch_limits_msTG = epoch_limits_msTG ; anapar.epoch_limits_msR = epoch_limits_msR; anapar.blinkTh = blinkTh; anapar.artifTh = artifTh; anapar.ARchecklim = ARchecklim; anapar.ARchans =ARchans;
    save([bigmatfolder allsubj{s} '_intp'],'erp','t','modir','coh','cond','respLR','RT','blocknum','filenames','anapar','artifact','erpr','validrlock','blink','artifact');
end
% spruce up this figure, which orders the maxmaxabs across trials and gives a sense of the proportion of trials that
% would make the cut, AND how that might potentially change if you were to adjust the threshold
figure(cumMaxAbsFig1); title('ordered maxmaxabs for each subject'); legend(allsubj(1:s)); xl=xlim; plot(xl,[1 1]*artifTh,'k--') % will need to zoom in!
ylabel('maxmaxabs'); xlabel('proportion of trials')
figure; set(gcf,'Position',[100 180 1000 400]); subplot(1,3,1), bar(propBlink); title('propBlink'),  subplot(1,3,2), bar(propArtifact); title('propArtifact'); subplot(1,3,3), bar(propKept); title('propKept'); 

% Much better! Proportions of artifact-free trials kept are higher now. Looking at individual histograms, there is one case, P09, where it
% seems that even on trials with no blink, electrodes 69 and 93 are disproportionately responsible for trial rejection, so go back and add
% those to ch2interp for that subject. It only really buys us a handful more trials.

% Now zooming in to the ordered maxmaxabs figure, it looks like our chosen artifact rejection threshold etc seem fine. One subject is going
% to have a bit more rejection, but the rest are <30% rejection.
% So now we go up and uncomment the last bit of that loop that re-saves the _intp files with the flags for artifacts and blinks, and we can
% use those in the analyses later.