% This script now goes through each subject and derives ERPs by averaging across all the "good" trials for each relevant condition
% of each subject, and plots the stimulus-locked and response-locked ERPs

clear all
allsubj = {'P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16'}; 

% the usual stuff:
bigmatfolder = 'bigmats/';
load chanlocsBioSemi128;
nchan = 128; % number of channels
next = 8;  % number of externals
fs=1024;
VEOGchans = [129 130]; 

load([bigmatfolder 'P01_intpCSD'],'anapar') % read in the analysis parameters that were used - be sure they were the same for all subjects!
epoch_limits_msTG = anapar.epoch_limits_msTG; % [-53 747];    % Target-locked epoch. We could use the same windows as were used for artifact rejection, or we might want different windows here - if so, don't forget that artifacts were only checked for in the window specified in IdentifyBadTrials!
tts = round(epoch_limits_msTG(1)/1000*fs):round(epoch_limits_msTG(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
tt = tts*1000/fs; % hence timebase in milliseconds, for plotting etc
epoch_limits_msR = anapar.epoch_limits_msR; % [-427 53];    % Response-locked epoch; again using integer number of SSVEP cycles (1000/18.75 = 53.33)
trs = round(epoch_limits_msR(1)/1000*fs):round(epoch_limits_msR(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
tr = trs*1000/fs; % hence timebase in milliseconds

clear avERP avERPr
for s=1:length(allsubj)
    
    subjID = allsubj{s};
    disp([subjID '...'])
    
    load([bigmatfolder subjID '_intpCSD'])

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
    
      
    % Get the average waveforms
    for c=1:2   % coherence
        % Select trials for averaging:
        trl = find(coh==c & ~blink & ~artifact & validrlock); % find the trial indices that had this particular coherence and  this evidence strength (or whatever YOU are looking for), and have no blinks or artifacts
        numtr(c,s) = length(trl); % record the number of good trials that went into each average
        avERP(:,:,c,s) = mean(erp(:,:,trl),3);
        avERPr(:,:,c,s) = mean(erpr(:,:,trl),3);
    end

    propKept(s) = length(find(~blink & ~artifact))/length(RT); 
end

%%

% select certain subjects to plot?
% sbj = 1:16; % usually you would plot all subjects like so
sbj = find(propKept>.7) % ... but you MAY want to exclude subjects who have a low trial count like this

% Important first of all to examine where activity shows up on the scalp during the trial
% Plot a series of topographies centered on these times:
TT = [0 160 320 480 640]; % you might be interested in more fine-grained steps if you have a sudden change at target onset, and want to see the series of evoked potential component caused by that, and assess their potential overlap and interference with the component we're mainly interested in here, the CPP. This is quite important to do. Here the taret onset was gradual so I'm spacing out the time points. 0 is good to include as a sanity check - the scalp should be pretty inactive.In discrete paradigms you have anticipatory processes already visible at t=0ms though.
winsize = 160; % window size in ms
figure; 
scale = [-1 1]*20; % sometimes you just take a punt at the scale (defines what value is deepest blue and what deepest red) and look at the result and then redefine to get it right. You should make sure the full range of colour is being filled, to best observe all the components happening
% we'll plot low coherence in top row, high coherence in bottom row
for i=1:length(TT)
    subplot(2,length(TT),i)
    trange = find(t>TT(i)-winsize/2 & t<TT(i)+winsize/2); % get the time indices covering the desired time range, centred on each point in TT
    topoplot(double(mean(mean(avERP(1:nchan,trange,1,sbj),4),2)),chanlocs,'electrodes','labels','colormap','jet','maplimits',scale) % some additional arguments that you can add: ,'maplimits',[-1 1]*4 sets the colorscale
    subplot(2,length(TT),length(TT)+i)
    topoplot(double(mean(mean(avERP(1:nchan,trange,2,sbj),4),2)),chanlocs,'electrodes','labels','colormap','jet','maplimits',scale) % some additional arguments that you can add: ,'maplimits',[-1 1]*4 sets the colorscale
    title([num2str(TT(i)) ' ms'])
end
% This shows a gradually emerging CPP among other things

% Same for response locked:
TTr = [-8:2:0]*1000/18.75; % stepping two ssvep cycles at a time
winsize = 2*1000/18.75; % window size in ms
figure; 
% we'll plot lowcoherence in top row, high coherence in bottom row
for i=1:length(TTr)
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
yl = [-2 28]; % y axis limits - consistent for Stimulus-locked and Response-locked
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

% Your figures for the first 16 subjects should look like RlockTopos_waveforms.png