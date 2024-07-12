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

% Let's say we also want to compute spectral anmplitude as a function of time in the Mu/beta bands (reflects motor preparation when measured
% over motor cortices). One way to do such time-frequency analysis is the short-Time Fourier Transform - taking FFTs in a sliding window
% across time. First define parameters of this analysis:
fftlen = round(fs/18.75*6); % Window of how many sample points? If there is an SSVEP involved, whether or not you are interested in analyzing it, it is good to have all power related to the SSVEP isolated in a single frequency bin. This happens when you choose a window length that is an integer number of SSVEP cycles. 
F = [0:fftlen-1]*fs/fftlen; % frequency scale, given window length (remember resolution = 1/window-duration)
ff = find((F>8 & F<18.7) | (F>18.8 & F<30)); % the indices of F that cover the spectral range/band of interest. Let's say we're interested in Mu and Beta bands combined. Note here I'm avoiding the SSVEP frequency (18.75hz in this example) 
Ts = [-300:50:1000]; % in msec, centered on what times do you want to measure spectral amplitude? i.e., where to center each consecutive window in time
Tr = [-450:50:0]; % for response-locked

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

    % Now in the rest of the code in the loop, we turn to computing and averaging Mu/Beta amplitude
    % compute short time fourier transform (STFT) for each single trial:
    STFT = []; % initialise
    % for the stimulus-locked STFT, we can compute ffts across all trials at once - that's why there is no 'for n=1:ntr' loop.
    for tt=1:length(Ts) % for each STFT timepoint (window centre) we'll compute the FFT for all trials at once
        [blah,samp] = min(abs(t-Ts(tt))); % find the sample point in the ERP epoch corresponding to the centre of the current FFT window
        spec = abs(fft(erp(:,samp-round(fftlen/2)+[1:fftlen],:),[],2))./(fftlen/2); % compute the magnitude of FFT (and scale it so that it's amplitude in same units as EEG
        % Save the result for each trial, just like the matrix 'erp'
        STFT(:,tt,:) = mean(spec(:,ff,:),2);
    end
    % response-locked (again processing each single trial):
    STFTr = nan(size(erp,1),length(Tr),ntr); % initalise. This time to NaN, because we want to retain same dimensions (a layer for each single trial) for indexing, but if R is such that we can't chop out a full STFT timecourse for a given trial, we mark that as NaN
    validrlockSTFT = zeros(1,ntr); % Again, make vector indicating whether valid for response locking (1) or not (0). The criteria for this are adjusted to account for windowing in STFT
    for n=1:ntr
        [blah,RTsamp] = min(abs(t-RT(n)*1000)); % find the sample point closest to the RT for this trial. Remember RT is in sec
        if RTsamp+Tr(1)*fs/1000-fftlen/2>0 & RTsamp+Tr(end)*fs/1000+fftlen/2<=length(t) & ~isnan(RT(n)) % Again, a trial is only validrlock=1 if the RT lies within the erp epoch you're using 
            validrlockSTFT(n)=1;
            for tt=1:length(Tr)
                spec = abs(fft(erp(:,RTsamp+round(Tr(tt)*fs/1000)-round(fftlen/2)+[1:fftlen],n),[],2))./(fftlen/2); % compute the magnitude of FFT (and scale it so that it's amplitude in same units as EEG
                STFTr(:,tt,n) = mean(spec(:,ff),2);
            end
        end
    end

    % Get the average Mu/Beta waveforms:
    % Here, note that motor lateralisation is critical, so you have to separate out the two alternatives, left and right, indexed by m. 
    % But note the trial designation of L/R can mean different things: whether the L or R button was PRESSED on a given trial, or whether L/R was correct on the
    % trial! (in decision-bias studies with prior cues, it could also mean whether the cue pointed L or R!)
    % You also have to decide whether you're averaging only correct trials (in which case the above distinction is moot), or all trials regardless of correctness.
    % My rules of thumb are that if you want to be able to measure differences in motor prep starting point across different blocked
    % regimes, you should compute Stimulus-locked Mu/Beta by taking all trials regardless of correctness, and if you want to see how buildup
    % is affected by different conditions, separate the trials by L/R motion direction (i.e. correct side), so you can separately plot the
    % motor prep signals contralateral vs ipsilateral to the correct side. 
    % However, for the R-locked, often the priority is to test whether there is a stereotyped level reached at response, and in that case
    % what you want is to separate trials by which button was pressed (whether correct or error). We'll do that here.
    for c=1:2   % condition
        for m=1:2 % direction of motion (or for R-locked, hand of response)
            % Select trials for averaging:
            trl = find(coh==c & modir==m & ~blink & ~artifact & validrlockSTFT);
            avMB(:,:,m,c,s) = mean(STFT(:,:,trl),3);

            % and response-locked: 
            trl = find(coh==c & respLR==m & ~blink & ~artifact & validrlockSTFT); % again, I am not selecting for correctness - you have to think about whether you want to do that
            avMBr(:,:,m,c,s) = mean(STFTr(:,:,trl),3);
        end
    end
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

% Your figures for the first 16 subjects should look like RlockTopos_waveforms.png. 

%% Now let's look at Mu/Beta 'MB'

% Let's first verify that there is lateralisation of MB, where the amplitude is lower contralateral to the button that was pressed. Plotting
% this topography also serves to highlight which electrodes might be best for measuring MB (although note there are some tasks like the
% continuous dots (2013 paper) and any of our delayed-response tasks, where there is precious little difference in MB amplitude contra/ipsi
% to responding hand, i.e. not much lateralisation):

trange = find(Tr>=-50,1); % pick a time just before response
figure;
topoplot(double(mean(mean(mean(avMBr(1:nchan,trange,1,:,sbj),5),4),2)-mean(mean(mean(avMBr(1:nchan,trange,2,:,sbj),5),4),2)),chanlocs,'electrodes','labels','colormap','jet');%,'maplimits',scale) % some additional arguments that you can add: ,'maplimits',[-1 1]*4 sets the colorscale
title('Left minus Right press')
colorbar
%% Mu/beta waveforms

ch = [3*32+19 32+22]; % select left/right channels - typically D19 and B22, and that's the case here
figure; hold on
for c=1:2
    plot(Ts,(mean(avMB(ch(2),:,1,c,sbj),5)+mean(avMB(ch(1),:,2,c,sbj),5))/2,'Color',colours(c,:),'LineWidth',2) % contralateral to correct side
    plot(Ts,(mean(avMB(ch(1),:,1,c,sbj),5)+mean(avMB(ch(2),:,2,c,sbj),5))/2,'--','Color',colours(c,:),'LineWidth',2) % ipsilateral (dashed)
end
set(gca,'Ydir','reverse') % we often turn the y axis upside down, to show increasing motor preparation (which is reflected in decreasing MB amplitude)

% Lateralisation - set it up so upwards means more preparation for the correct alternative
figure; hold on
for c=1:2
    plot(Ts,(mean(avMB(ch(1),:,1,c,sbj),5)+mean(avMB(ch(2),:,2,c,sbj),5))/2-(mean(avMB(ch(2),:,1,c,sbj),5)+mean(avMB(ch(1),:,2,c,sbj),5))/2,'Color',colours(c,:),'LineWidth',2)
end

% Response-locked:
figure; hold on
for c=1:2
    plot(Tr,(mean(avMBr(ch(2),:,1,c,sbj),5)+mean(avMBr(ch(1),:,2,c,sbj),5))/2,'Color',colours(c,:),'LineWidth',2) % contralateral to SIDE OF RESPONSE
    plot(Tr,(mean(avMBr(ch(1),:,1,c,sbj),5)+mean(avMBr(ch(2),:,2,c,sbj),5))/2,'--','Color',colours(c,:),'LineWidth',2) % ipsilateral (dashed)
end
set(gca,'Ydir','reverse') % we often turn the y axis upside down, to show increasing motor preparation (which is reflected in decreasing MB amplitude)

% notice that the traces coalesce at the time of response, as we've reported lots of times.