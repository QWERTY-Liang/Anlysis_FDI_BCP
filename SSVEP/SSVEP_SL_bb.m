%%  6.1 SSVEP here no baseline
% Author: Liang Tong
% Date: 24/7/2024 

%to update list
%1.function all the code to make main script simplier
%% Toolbox requirements: 
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
 EEG = pop_loadset([exp.finalpath 'nobaseline_SL_bTLalldata.set']);
 load TL_AllBehaviour_SL_nobaseline.mat 
 AllBehaviour_SL_bb=AllBehaviour_SL_b;
 clear AllBehaviour_SL_b
%% variable prepare
fs=EEG.srate;%512
t=EEG.times;
erp=EEG.data;
fv = 21.5; % flicker frequency in Hz - this is for EACH orientation. the overall on-off flicker rate of the stimulus will be double this
nchan=128;
chanlocs=EEG.chanlocs;
clear EEG

%% Now we'll do a spectral analysis on each single trial to pull out the timecourse of the 18.75Hz component

% We can use any spectral decomposition method - wavelets, short-time Fourier transform (STFT), etc. 
% Here we'll use something equivalent to an STFT as we typically have in past studies. An STFT computes an FFT in a
% certain short-ish time window, then shifts the window along in time a little and does it again, and again, etc, and
% then you can extract a certain frequency or frequency band from each fft and you have a time-course of its spectral
% amplitude.
% FFTs are convolutions of a signal with a complex sinusoid of a series of frequencies. We only care about the flicker frequency. 
% So to save space let's make a complex sinusoid at just 18.75 Hz to convolve with our EEG. We will pick a time window
% that has an integer number of flicker cycles in it, for best accuracy:
%fftlen = 164; % 273;  % in sample points. 164 samples is 320 ms, almost exactly 6 cycles of the flicker.
fftlen = round(fs/21.5*6);
% Make the complex sinusoid. It's rather like a wavelet but is not tapered like normal wavelets are:
clear exp
wvlt = exp(i*2*pi*fv*[1:fftlen]/fs);

% now convolve this with all erps at all channels:
clear sv
for e=1:nchan
    disp(e);% show the process
    for n=1:size(erp,3)
        sv(e,:,n) = conv(erp(e,:,n),wvlt,'same');
    end
    
end
% (takes a while)

% After doing this, sv gives us a complex 18.75 Hz signal for each channel for each single trial.
% The absolute phase of sv depends on things like transmission delays to visual cortex as well as the choice of phase of the complex 'wavelet', and it doesn't really matter. 
% What matters is the relative phase among the signals. 

%% Now let's plot some averages across trials in various conditions
% In the task code, we pseudo-randomised which grating came first in the flicker.
% The indicator vector firstTilt says whether the first grating presented was a left-tilted grating (1) or a right-tilted one (2).
% So firstTilt=1 means left-off-right-off-left-off-ri...
% &  firstTilt=2 means right-off-left-off-right-off-lef...
% From this we can compute whether the first grating was the brighter one, regardless of tilt:

 par.PTonecycleL = [1,0,0,0;0,0,1,0];  % pulse train for one cycle for left-tilted, flipped also
 par.PTonecycleR = [0,0,1,0;1,0,0,0];  % pulse train for one cycle for right-tilted
firstTilt=AllBehaviour_SL_bb(:,10);% left=1 right=2
brighterLR= AllBehaviour_SL_bb(:,9);% left right=2
durationCoh=AllBehaviour_SL_bb(:,3);% 0.14 and 0.07
respLR=AllBehaviour_SL_bb(:,8);
muscle=AllBehaviour_SL_bb(:,4);

% 
 brightfirst = (firstTilt==brighterLR); % 0 is first dark; 1 is first bright

% In this particular experiment we have 5 duration/coherence combinations, and 2 possible 'directions' (which is
% brighter, left or right-tilted), so 10 conditions altogether.
% durationCoh = 1 means low contrast difference and duration 0.2 sec
% durationCoh = 2 means low contrast difference and duration 0.4 sec
% durationCoh = 3 means low contrast difference and duration 0.8 sec
% durationCoh = 4 means low contrast difference and duration 1.6 sec
% durationCoh = 5 means high contrast difference and duration 1.6 sec
% Critically, the duration of the STIMULUS is always 1.6 sec beyond the evidence onset, so in durationCoh=2 for example, the
% grating contrasts step back to 50-50 after 400 ms.

% Get the average waveforms for each condition regardless of outcome:
%% High vs low
clear SV
contrast=[0.07 0.14];
for c=1:2  % durationCoh condition
    for lr=1:2 % which grating was brighter?
        for bf=0:1 % brightfirst
            trl = find(durationCoh==contrast(c) & brightfirst==bf & brighterLR==lr & AllBehaviour_SL_bb(:,5)~=3 & AllBehaviour_SL_bb(:,5)~=5 );% & respLR~=brighterLR);% & ~blink & ~artifact);
           
            numtr(c,lr,bf+1) = length(trl); % good to record the number of trials that went into each average
            SV(:,:,c,lr,bf+1) = mean(sv(:,:,trl),3);
        end
    end
end
% note there's no artifact rejection - that should be added as you prefer!



%%
% The most obvious thing we should expect of course is that the
% signal has opposite phase for brightfirst==0 relative to brightfirst==1. 
% Let's have a look at electrode 23 (Oz) - ssvep for a stimulus centred on the fovea should usually show up around there:
% ch=23;
% figure; hold on
% col = jet; % choose some colours for the traces
% colors = col([1:15:64],:); 
% 
% % Now SV are complex exponentials which are effectively two-dimensional (real and imaginary part), 
% % so to look at it as a simple sinusoidal signal over time, we plot just the real part (or the imaginary part, again all that matters is RELATIVE phase betwen conditions)
% 
% for c=1:2
%     % We'll just average across brightLR =1/2 presuming that's symmetric
%     plot(t,real(mean(SV(ch,:,c,:,1),4)),'Color',colors(c,:),'LineWidth',1) % Dim first
% end
% 
% for c=1:2
%     plot(t,real(mean(SV(ch,:,c,:,2),4)),'Color',colors(c,:),'LineWidth',2) % Bright first
% end
% legend(num2str([1:2]'))
ch = 23;
figure;

% Choose colors for conditions
colors = [0 0 1; 1 0 0]; % Blue for condition 1, Red for condition 2

% Define the time points and labels for xline
xlines = [-1300, -1200, -600, 0, 800, 1500, 2000];
xlabels = {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s', 'DDL-2s'};

% Dim first subplot
subplot(2, 1, 1); hold on;
title('Dim First');

for c = 1:2
    plot(t, real(mean(SV(ch, :, c, :, 1), 4)), 'Color', colors(c, :), 'LineWidth', 1);
end

% Add vertical lines with labels
for i = 1:length(xlines)
    xline(xlines(i), '--r', xlabels{i});
end

% Bright first subplot
subplot(2, 1, 2); hold on;
title('Bright First');

for c = 1:2
    plot(t, real(mean(SV(ch, :, c, :, 2), 4)), 'Color', colors(c, :), 'LineWidth', 1);
end

% Add vertical lines with labels
for i = 1:length(xlines)
    xline(xlines(i), '--r', xlabels{i});
end

% Add legends to the subplots
subplot(2, 1, 1);
legend('Condition 1', 'Condition 2');

subplot(2, 1, 2);
legend('Condition 1', 'Condition 2');
% Note bright first and dim first are opposite phase as expected. The time it takes for amplitude to ramp up is about
% 320 ms, the fft window length, which makes sense. 
%% Since there is that expected tendency for brightfirst and dimfirst to be opposite phase, we can just subtract them:
figure; hold on
for c=1:2
    plot(t,real(mean(diff(SV(ch,:,c,:,:),[],5),4)),'Color',colors(c,:)) % brightfirst minus dimfirst
    %plot(t,real(mean(mean(SV(ch,:,c,:,:),5),4)),'Color',colors(c,:)) % brightfirst minus dimfirst
end
% Define the time points and labels for xline
xlines = [-1300, -1200, -600, 0, 800, 1500, 2000];
xlabels = {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s', 'DDL-2s'};
% Add vertical lines with labels
for i = 1:length(xlines)
    xline(xlines(i), '--r', xlabels{i});
end
legend('Condition 1', 'Condition 2');
% You can see from the high contrast difference condition that the ramp-up reaches its midpoint at about 100 ms and
% ramp-down midpoint at about 1700 ms, together suggesting there's a circa 100-ms delay from changes on the screen to
% changes in the Steady-state VEP amplitude. The shortest duration condition 1 has 200 ms of contrast difference so the
% SV amplitude should reach its highest at the midpoint plus 100 ms delay = 200 ms, which it does, and as expected, the
% 400-ms duration condition is the next to reduce in amplitude, then the 800 ms, but there are some very interesting
% dynamics beyond that, due to the subject actively analysing those grating contrasts to compare them.
% So next let's try to recapitulate these plots but taking a view where positive values mean towards the "correct" phase
% (i.e. driven by the grating that is actually higher in contrast on a given trial) and negative mean towards the
% opposite, "wrong" phase favouring the dimmer grating.
%%
% The way we'll do it is by choosing one reference timepoint where we most confidently believe the phase is driven by
% the true contrast difference, note the phase at that time, taking it as representative of the correct phase, and
% extrapolate back and forward from there.
    
% Seems most sensible to take a timepoint around the time the ssvep initially levels off (which is dependent on the fftlen of course)
% - around 250 ms. In computing the reference phase I'll just include conditions 2-4 because condition 1 is already on the way back down by 300 ms, 
% and condition 5 potentially has a very different phase due to having very different contrast difference (remember visual cortical response latencies scale inversely with contrast)
% Note that in using conditions 2-4 as reference phase, we'll likely be underestimating the true amplitude of the
% high-contrast condition so if we really cared about optimally extracting the high-contrast condition we'd compute and 
% use a different reference phase for that one condition. But I'm not bothering with that here.
% 
% reftime = find(t>500,1); % reference time
% refphase = angle(mean(mean(diff(SV(ch,reftime,1:2,:,:),[],5),4),3)); % reference phase 
% %refphase = angle(mean(mean(mean(SV(ch,reftime,1,:,:),5),4),3)); % reference phase 
% % now make a reference complex exponential wave that extrapolates from that reference time/phase:
% refwave = exp(i*((t'-t(reftime))/1000*fv*2*pi)+refphase);
% 
% plot(t,real(refwave)*500,'k') % can plot it in the same way as the averages above to make sure it lines up at the
% % reference time point.
% 
% % Now use this reference rotate all complex exponential signals at all time points for all single trials so that the values represent the 
% % degree to which the differential ssvep leans towards the correct, higher-contrast-driven phase, by extrapolating from this reference phase:
% % We have to rotate the dimfirst trials in the opposite direction so we'll do the two separately:
% clear sv1
% sv1(:,find(brightfirst==1)) = squeeze(sv(ch,:,find(brightfirst==1)))./repmat(refwave,[1,length(find(brightfirst==1))]);
% sv1(:,find(brightfirst==0)) = squeeze(sv(ch,:,find(brightfirst==0)))./repmat(refwave,[1,length(find(brightfirst==0))]);
% 
% % For reasons best explained in secondary school trigonometry, dividing by a unit-amplitude complex exponential has the
% % effect of rotating the angle backwards so that for all time points it should tend to lie on the positive real axis.
% % All we have to do below then is take the real part of the signal, and effectively, what we're seeing is the ssvep
% % amplitude "projected" onto the "correct" phase.
% 
% %% Now we can divide up and plot averages across trials in whatever way we want 
% % because for each single trial, positive means toward the higher-contrast grating and negative means away
% clear SV1
% for c=1:2   % durationCoh condition
%     SV1(:,c) = mean(sv1(:,find(durationCoh==contrast(c) & brightfirst==1)),2);
% end
% 
% figure; hold on
% for c=1:2
%     plot(t,real(SV1(:,c)),'Color',colors(c,:)) 
% end
% legend(num2str([1:2]'))
% 
% % what we can see here is that when the evidence switches off and goes back to equal-contrast, there is a period of time
% % where the grating that WAS dimmer wins out (negative sv1 values), probably due to selective adaptation, and then
% % recovers to around zero. One of our hypotheses for this dataset was that there might be a confirmation bias expressed
% % in the sensory evidence signal, where, because some evidence has been gathered in favour of the grating that was
% % brighter for some time, that grating is boosted even after it drops back down to the baseline contrast. Not happening
% % for this subject anyway!
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% New (simpler) approach:
% 
% 
% % turn all brightfirst = 0 upside down to bring into phase alignment with brightfirst = 1, for convenience:
% sv(:,:,find(brightfirst==0)) = -sv(:,:,find(brightfirst==0));
% % Now instead of extrapolating phase from one time point, we will use the 1.6 sec duration low contrast condition as the
% % reference:
% clear SVref4
% SVref4 = mean(sv(:,:,find(durationCoh==0.07)),3);
% refphase = angle(SVref4);
% % now rotate all trials so that if they align with this condition 4, they will have average phase = 0 (i.e. rotate around to align with positive real axis in complex plane)
% sv1 = sv./repmat(exp(i*refphase),[1,1,size(sv,3)]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. High vs Low(correct vs wrong)
% Now make some averages:

clear SV1 SV1corr SV1err
for c=1:2  % durationCoh condition
    SV1(:,:,c) = mean(sv1(:,:,find(durationCoh==contrast(c))),3);
    SV1corr(:,:,c) = mean(sv1(:,:,find(durationCoh==contrast(c)& respLR==brighterLR)),3);
    SV1err(:,:,c) = mean(sv1(:,:,find(durationCoh==contrast(c) & respLR==3-brighterLR)),3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. FDI vs BCP
clear SV1 SV1corr SV1err
for c=1:2  % durationCoh condition
    SV1(:,:,c) = mean(sv1(:,:,find(muscle == c)),3);
    SV1corr(:,:,c) = mean(sv1(:,:,find(muscle == c & respLR==brighterLR)),3);
    SV1err(:,:,c) = mean(sv1(:,:,find(muscle == c & respLR==3-brighterLR)),3);
end

%% Topographies
% Now a topography to pick a 'best' electrode: we can do this based on SIGNAL TO NOISE RATIO of condition 4 (because we 
% can average across the longest interval for robustness), taking the baseline SSVEP as the noise level:
topodat_SN = mean(abs(SVref4(:,find(t>300 & t<1500))),2) ./ mean(abs(SVref4(:,find(t>-500 & t<0))),2);
figure; topoplot(double(topodat_SN),chanlocs,'electrodes','labels'); colorbar

% electrode 9 has best SNR for this subject
[~,I] = sort(topodat_SN,'descend');
Bestelectrode_SN = I(1);
% Topographies of correct-error difference:
trange = find(t>300 & t<800); % taking last half second because that seems to be when many subjects separate their correct and error
figure;set(gcf,'Position',[50 80 1800 500]);
for c=1:2
    subplot(1,2,c)
    topodat_CorrErr = real(mean(SV1corr(:,trange,c),2)) - real(mean(SV1err(:,trange,c),2));
    [~,I] = sort(topodat_CorrErr,'descend');
    Bestelectrode_CorrErr(c) = I(1);
    topoplot(double(topodat_CorrErr),chanlocs,'electrodes','labels','colormap','jet'); 
end

% so the electrode that has strongest relationship with choice is electrode 22! 
% Now electrode 9 will have best 'neurometric performance' in that it is most strongly related to true contrast
% condition, but electrode 22 has the strongest relationship with choice, so is probably the one that makes more sense
% to look at impact on the CPP. But we should try both!
%% plot traces at these two electrodes to have a look
col = jet; % choose some colours for the traces
%colors = col([1:33:180],:);
colors = col([1:50:256],:);

for e = [Bestelectrode_SN Bestelectrode_CorrErr(1)]%[9 22] 
    figure; set(gcf,'Position',[50 80 1200 500]);
    %SSVEP
    subplot(1,2,1); hold on
    for c=1:2, plot(t,real(SV1(e,:,c)),'Color',colors(c,:),'LineWidth',1), end
    plot([t(1) t(end)],[0 0],'Color',[.6 .6 .6]); legend(num2str([1:2]'))
    title(['elec' num2str(e) ': SSVEP ' ])
    % correct vs error:
    %set(gca,'FontSize',14)

    subplot(1,2,2); hold on
    for c=1:2
        plot(t,real(SV1corr(e,:,c)),'Color',colors(c,:),'LineWidth',1)
        plot(t,real(SV1err(e,:,c)),'--','Color',colors(c,:),'LineWidth',1)
    end
    plot([t(1) t(end)],[0 0],'Color',[.6 .6 .6]); 
    title(['elec' num2str(e) ': SSVEP CorrErr '])
    
end


% It is a good idea to check for differences in phase and amplitude together using a polar plot where we can plot the
% conditions as full complex numbers. Note that for electrode 9 you can see that there is a somewhat consistent phase
% lag for error relative to correct
for e=[Bestelectrode_SN Bestelectrode_CorrErr(1)] %[9 22]
    figure; polar(max(SV1corr(:))); hold on
    trange = find(t>1000 & t<1500);
    for c=1:2
        plot(mean(SV1corr(e,trange,c),2),'o','Color',colors(c,:),'LineWidth',2)
        plot(mean(SV1err(e,trange,c),2),'d','Color',colors(c,:),'LineWidth',2)
    end
    title(['elec ' num2str(e) ' Copmlex Corr/Err'])
    
end
