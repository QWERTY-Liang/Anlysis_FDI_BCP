% This script now goes through each subject and derives average ERPs for relevant conditions of each subject, 
% and plots the stimulus-locked and response-locked ERPs

clear all
allsubj = {'P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'RL01', }; 

bigmatfolder = 'bigmats/';
load chanlocsBioSemi128;
nchan = 128; % number of channels
next = 8;  % number of externals

fs=1024;

% Before sorting trials and averaging by condition, we need to figure out what analysis time window and artifact
% rejection (AR) criteria will give us the best data quality across as many subjects as possible.
% AR criteria by our very simply threshold method consist simply of a time window in which to check for
% artifacts, and a threshold level that marks a trial as having an artifact if exceeded within that window.
% The first thing to consider is blink behaviour because it's the biggest artifact source. If you instructed
% the subjects well, they will have waited until the stimulus has fully ended or until a good 0.5-1 sec after
% responding before they get their blink in, but sometimes they blink very quickly and if you set too wide an
% AR window in this case, loads of trials will be rejected even if the time between stimulus and response is
% mostly clean. When there are several such subjects in a dataset, it can be best to artifact reject just
% based on a relatively narrow window from stimulus onset and back from the response - this is what I arrive
% at below in fact. It is far better to instrct your subjects well than to have to scrape around for trials by
% adjusting windows like this!

VEOGchans = [129 130]; %  CHECK THAT THESE ARE INDEED YOUR VEOG CHANNELS [upr lwr]
%%
% We will now 'draft' some artifact checking settings and see how many trials get rejected, while also looking
% at the blink behaviour from the VEOG channels
blinkTh = 100; % threshold for detecting a blink in VEOG upper-lower difference signal
artifTh = 700; % this is for SCALP channels - it's high because it is in CSD units, which take bigger values
ARchans = [1:128]; % Artifact rejection channels - which electrodes should we reject based on?
ARchecklim = [-100 1760]; % in msec, the limits of the time window relative to the event in which you will check for an artifact

blinkTms_all = []; blinkTms_allR = []; % (initialisation) these will contain the times relative to stimulus and to response ('R') at which blinks occurred on the trials
for s=1:length(allsubj)
    
    subjID = allsubj{s};
    disp([subjID '...'])
    
    load([bigmatfolder subjID '_intp'])

    % First analyse blinks. Let's assume most of the time window has eyes-open so the median value will
    % represent that state, to be used as a baseline:
    VEOG = squeeze(erp(VEOGchans(1),:,:)-erp(VEOGchans(2),:,:)); % Upper elec minus lower elec; so now it's time x trials
    VEOG = VEOG - repmat(median(VEOG),[size(VEOG,1),1]); % "subtract the value estimated as the eyes open state, so blinks are seen as deviations from this
    % now plot VEOG to investigate blinking behaviour:
    figure; plot(t,VEOG); title(['Subject ' allsubj{s}]) % plot the VEOG waveforms for all trials
    % For each trial get the time at which the blink threshold is exceeded, to look at the distribution of
    % times subjects tend to blink:
    blinkTms=[]; 
    for n=1:size(VEOG,2), 
        blnk=find(VEOG(:,n)>blinkTh,1); 
        if isempty(blnk), blinkTms(n)=nan; else, blinkTms(n)=t(blnk); end; 
    end; 
    figure; subplot(1,2,1); hist(blinkTms,[t(1):100:t(end)]); title(['Subj ' allsubj{s} 'S-lock']); subplot(1,2,2); hist(blinkTms-RT*1000,[-1000:100:1000]); title('R-lock')
    blinkTms_all = [blinkTms_all blinkTms]; % and append to a big vector of blink times to check across all subjects
    blinkTms_allR = [blinkTms_allR blinkTms-RT*1000];
    
    % the following will count up all of the trials with artifacts detected, whether blinks or other noise on
    % the cap channels:
    blink = max(abs(VEOG)) > blinkTh; % this will give a logical 1 if this threshold is exceeded on a certain trial, 0 otherwise
    % and then artifacts on any of the channels specified in ARchans above
    artcheckwin = find(t>=ARchecklim(1) & t<=ARchecklim(2)); % indices of the timepoints for the artifact check window whose limits (min and max) are defined above
    art = squeeze(max(abs(erp(ARchans,artcheckwin,:)),[],2)-artifTh > 0);  % (:,artcheckwin,n) % art is now elec x trial, 1 if any timepoint exceeds artif thresh, 0 if fine
    numart = sum(art,2); % How many trials had an artifact for each ELECTRODE    
    artifact = sum(art,1); % How many electrodes had an artifact on each TRIAL
    figure; bar(numart); title(['Subject ' allsubj{s}])
    
    propKept(s) = length(find(~blink & ~artifact))/length(RT); % proportion of trials that would be kept if these were the artifact/blink detection criteria
end
propKept
%%
% First of all, the proportion of trials not rejected based on either a blink or a large amplitude on any cap
% electrode in the long artifact check window ranges from excellent (95%) to terrible (3%). Distribution:
figure; hist(propKept) % - 8 subjects with <50% trials retained. It's possible there are that many bad subjects 
% but it's worth trying to tighten up the artifact check window

% As you can see from the VEOG plots, some subjects are good, some are terrible blinkers. For many subjects you can see that they
% rarely blink during the short time between evidence onset and 0.5-1 sec but relative to the response, they
% blink very soon afterwards. Let's check the histograms of blink times across all subjects (pooled):
figure; subplot(1,2,1), hist(blinkTms_all,[t(1):100:t(end)]); title('S-lock'); subplot(1,2,2); hist(blinkTms_allR,[-2500:100:1500]); xlim([-2500 1500]); title('R-lock')
    
% This is quite terrible - subjects  almost look like they are trying to time a blink on their response. With
% better instructed subjects you could choose a wide range of alternative artifact check strategies and still
% get reliable results but with this dataset, we will clearly have to carefully tailor the settings so as to
% squeeze as many trials as possible out of each subject...
% The goal is always to retain as many trials as possible but to avoid letting in noisy trials, to get the most reliable average
%%
% This is our plan:
% We will choose target-locked and response-locked windows that are just wide enough to capture the dynamics 
% of the decision process just after target onset and just before response, and only artifact-reject based on data in those narrow windows

epoch_limits_msTG = [-53 747];    % Target-locked epoch; again using integer number of SSVEP cycles (1000/18.75 = 53.33)
tts = round(epoch_limits_msTG(1)/1000*fs):round(epoch_limits_msTG(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
tt = tts*1000/fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds
epoch_limits_msR = [-427 53];    % Response-locked epoch; again using integer number of SSVEP cycles (1000/18.75 = 53.33)
trs = round(epoch_limits_msR(1)/1000*fs):round(epoch_limits_msR(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
tr = trs*1000/fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds
% Note that in Response-locked waveforms time '0' is the time of response, and weusually don't care so much about 
% activity after that time. Our strategy is to cut out the shorter T-locked and R-locked epochs from EACH of the 
% single stimulus-locked epochs. We have to do it on a trial by trial basis (inside the loop below) because RT is different from trial to trial!

clear avERP avERPr
for s=1:length(allsubj)
    
    subjID = allsubj{s};
    disp([subjID '...'])
    
    load([bigmatfolder subjID '_intp'])

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
    art = squeeze(max(abs([erp(ARchans,artcheckwinTG,:) erpr(ARchans,:,:)]),[],2)-artifTh > 0);  % art is elec x trial, 1 if any timepoint exceeds artif thresh, 0 if fine. Again, concatenated t-lock and R-lock so will detect artifact in EITHER
    numart = sum(art,2); % How many trials had an artifact for each ELECTRODE    
    artifact = sum(art,1); % How many electrodes had an artifact on each TRIAL
    figure; bar(numart); title(['Subject ' allsubj{s}])
      
    % Get the average waveforms
    for c=1:2   % coherence
        % Select trials for averaging:
        trl = find(coh==c & ~blink & ~artifact & validrlock); % find the trial indices that had this particular coherence and  this evidence strength (or whatever YOU are looking for), and have no blinks or artifacts
        numtr(c,s) = length(trl); % record the number of good trials that went into each average
        avERP(:,:,c,s) = mean(erp(:,:,trl),3);
        avERPr(:,:,c,s) = mean(erpr(:,:,trl),3);
    end

    propKept(s) = length(find(~blink & ~artifact))/length(RT); % proportion of trials that would be kept if these were the artifact/blink detection criteria
end
propKept % better? Quite. There are three subjects who have >30% trials rejected that should probably be excluded

%%

% select certain subjects to plot?
% sbj = 1:24; % all?
sbj = find(propKept>.7)

% Important first of all to examine where activity shows up on the scalp during the decision
% Plot a series of topographies centered on these times:
TT = [0 160 320 480 640];
winsize = 160; % window size in ms
figure; 
% we'll plot lowcoherence in top row, high coherence in bottom row
for i=1:length(TT)
    scale = [-16 16];
    subplot(2,length(TT),i)
    trange = find(t>TT(i)-winsize/2 & t<TT(i)+winsize/2); 
    topoplot(double(mean(mean(avERP(1:nchan,trange,1,sbj),4),2)),chanlocs,'electrodes','labels','colormap','jet','maplimits',scale) % some additional arguments that you can add: ,'maplimits',[-1 1]*4 sets the colorscale
    subplot(2,length(TT),length(TT)+i)
    topoplot(double(mean(mean(avERP(1:nchan,trange,2,sbj),4),2)),chanlocs,'electrodes','labels','colormap','jet','maplimits',scale) % some additional arguments that you can add: ,'maplimits',[-1 1]*4 sets the colorscale
    title([num2str(TT(i)) ' ms'])
end

% Same for response locked:
TTr = [-8:2:0]*1000/18.75; % stepping two ssvep cycles at a time
winsize = 2*1000/18.75; % window size in ms
figure; 
% we'll plot lowcoherence in top row, high coherence in bottom row
for i=1:length(TTr)
    scale = [-16 16];
    subplot(2,length(TTr),i)
    trange = find(tr>TTr(i)-winsize/2 & tr<TTr(i)+winsize/2); 
    topoplot(double(mean(mean(avERPr(1:nchan,trange,1,sbj),4),2)),chanlocs,'electrodes','labels','colormap','jet','maplimits',scale) % some additional arguments that you can add: ,'maplimits',[-1 1]*4 sets the colorscale
    subplot(2,length(TTr),length(TTr)+i)
    topoplot(double(mean(mean(avERPr(1:nchan,trange,2,sbj),4),2)),chanlocs,'electrodes','labels','colormap','jet','maplimits',scale) % some additional arguments that you can add: ,'maplimits',[-1 1]*4 sets the colorscale
    title([num2str(TTr(i)) ' ms'])
end

% Looking at topographies there is an interesting array of foci on the scalp, but the one that rises to the
% highest amplitude is focused on electrode A19 - we identify this as the CPP and plot waveforms below. 
% What are the other foci? Wel the posterior bilateral occipital negativity we often see during visual tasks
% and could relate to deploying resources to the visual system to scrutinise the target (very speculative);
% the frontal foci may relate to the CNV (though that's usually midline) and/or motor preparation

%% Plot waveforms against time:
ch = [19]; % choose channel to plot based on topographies above

% first create figure and set plot dimensions etc:
set(0,'DefaultLegendAutoUpdate','off')
figure; set(gcf,'Position',[80 80 700 400]); set(gcf,'DefaultLineLineWidth',2);
yl = [-2 26]; % y axis limits - consistent for Stimulus-locked and Response-locked
xlS = epoch_limits_msTG; % x axis limits for stimulus-locked
xlR = epoch_limits_msR;  % x axis limits for response-locked

fY=.14; % Position of bottom of plots within figure (these are all proportions of full figure size)
fW_S=.5; % plot width for stimulus-locked. response-locked will be scaled according to x-axis limits then, so that 100 ms spans the same physical distance in both panels
fH=.82; % plot height - should be same for both panels
FS = 14; % fontsize

% make both panels:
ax1 = axes; set(ax1,'Position',[.09 fY fW_S fH]); hold on
ax2 = axes; set(ax2,'Position',[.64 fY fW_S*diff(xlR)/diff(xlS) fH]); hold on

colours = [[120 120 255];[0 0 150]]/255; % light blue for low-coh, dark blue for high coh
% Now plot:
for c = 1:2
    axes(ax1); xlim(xlS); ylim(yl)
    plot(t,mean(mean(avERP(ch,:,c,sbj),4),1), 'color',colours(c,:), 'Linewidth', 2)
    axes(ax2); xlim(xlR); ylim(yl)
    plot(tr,mean(mean(avERPr(ch,:,c,sbj),4),1), 'color',colours(c,:), 'Linewidth', 2)
end
% Aesthetics:
axes(ax1);
plot([0 0],yl,'k','LineWidth',3)
xlabel('Time (msec)');
ylabel('CPP (\muV/m^{2})');
legend('40%','80%')
set(gca,'FontSize',14)
axes(ax2);
plot([0 0],yl,'k','LineWidth',3)
set(gca,'FontSize',14)
set(gca,'YTickLabel','');
set(gca,'YColor',[1 1 1]);