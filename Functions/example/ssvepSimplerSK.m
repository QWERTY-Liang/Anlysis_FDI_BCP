clc; clear all; close all;
fv = 18.75; % flicker frequency in Hz - this is for EACH orientation. the overall on-off flicker rate of the stimulus will be double this

% I forgot to tell you to save a variable called "firstTilt" for all single trials. I added it to the 
% matfile for S9 by using the simple code:
firstTilt = [];
for b=1:12 % 12 blocks , S1 10 blocks (1:10), S10(11 blocks, 1:11) S13(11, 2:12) S15(12, 2:13)
    load(['./S9_' num2str(b) 'VD'])
    firstTilt = [firstTilt par.firstTilt];
end
% save firstTiltS9 firstTilt
%%

subjID = 'S9'; % subject
load(['./' subjID '_intp'])
load(['./' subjID '_cpp'])

fs=512; % EEG sample rate

load chanlocsBioSemi128;
nchan = 128; % number of channels
next = 8;  % number of externals

%% Now we'll do a spectral analysis on each single trial to pull out the timecourse of the 18.75Hz component

% We can use any spectral decomposition method - wavelets, short-time Fourier transform (STFT), etc. 
% Here we'll use something equivalent to an STFT as we typically have in past studies. An STFT computes an FFT in a
% certain short-ish time window, then shifts the window along in time a little and does it again, and again, etc, and
% then you can extract a certain frequency or frequency band from each fft and you have a time-course of its spectral
% amplitude.
% FFTs are convolutions of a signal with a complex sinusoid of a series of frequencies. We only care about the flicker frequency. 
% So to save space let's make a complex sinusoid at just 18.75 Hz to convolve with our EEG. We will pick a time window
% that has an integer number of flicker cycles in it, for best accuracy:
fftlen = 164; % 273;  % in sample points. 164 samples is 320 ms, almost exactly 6 cycles of the flicker.
% Make the complex sinusoid. It's rather like a wavelet but is not tapered like normal wavelets are:
wvlt = exp(i*2*pi*fv*[1:fftlen]/fs);

% now convolve this with all erps at all channels:
% tic
clear sv
for e=1:nchan
    for n=1:size(erp,3)
        sv(e,:,n) = conv(erp(e,:,n),wvlt,'same');
    end
end

%% New (simpler) approach:

brightfirst = (firstTilt==brighterLR); 
% turn all brightfirst = 0 upside down to bring into phase alignment with brightfirst = 1, for convenience:
sv(:,:,find(brightfirst==0)) = -sv(:,:,find(brightfirst==0));
% Now instead of extrapolating phase from one time point, we will use the 1.6 sec duration low contrast condition as the
% reference:
clear SVref4
SVref4 = mean(sv(:,:,find(durationCoh==4)),3);
refphase = angle(SVref4);
% now rotate all trials so that if they align with this condition 4, they will have average phase = 0 (i.e. rotate around to align with positive real axis in complex plane)
sv1 = sv./repmat(exp(i*refphase),[1,1,size(sv,3)]);
%%
% Now make some averages:
clear SV1 SV1corr SV1err
for c=1:4   % durationCoh condition
    SV1(:,:,c) = mean(sv1(:,:,find(durationCoh==c)),3);
    SV1corr(:,:,c) = mean(sv1(:,:,find(durationCoh==c & respLR==brighterLR)),3);
    SV1err(:,:,c) = mean(sv1(:,:,find(durationCoh==c & respLR==3-brighterLR)),3);
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
trange = find(t>1000 & t<1500); % taking last half second because that seems to be when many subjects separate their correct and error
figure;set(gcf,'Position',[50 80 1800 500]);
for c=1:4
    subplot(1,4,c)
    topodat_CorrErr = real(mean(SV1corr(:,trange,c),2)) - real(mean(SV1err(:,trange,c),2));
    [~,I] = sort(topodat_CorrErr,'descend');
    Bestelectrode_CorrErr(c) = I(1);
    topoplot(double(topodat_CorrErr),chanlocs,'electrodes','labels'); 
end

% so the electrode that has strongest relationship with choice is electrode 22! 
% Now electrode 9 will have best 'neurometric performance' in that it is most strongly related to true contrast
% condition, but electrode 22 has the strongest relationship with choice, so is probably the one that makes more sense
% to look at impact on the CPP. But we should try both!
%% plot traces at these two electrodes to have a look
col = jet; % choose some colours for the traces
%colors = col([1:33:180],:);
colors = col([1:50:256],:);

for e = [Bestelectrode_SN Bestelectrode_CorrErr(4)]%[9 22] 
    figure; set(gcf,'Position',[50 80 1200 500]);
    %SSVEP
    subplot(1,2,1); hold on
    for c=1:4, plot(t,real(SV1(e,:,c)),'Color',colors(c,:),'LineWidth',2), end
    plot([t(1) t(end)],[0 0],'Color',[.6 .6 .6]); legend(num2str([1:4]'))
    title(['elec' num2str(e) ': SSVEP ' ])
    % correct vs error:
    set(gca,'FontSize',14)

    subplot(1,2,2); hold on
    for c=1:4
        plot(t,real(SV1corr(e,:,c)),'Color',colors(c,:),'LineWidth',2)
        plot(t,real(SV1err(e,:,c)),'--','Color',colors(c,:),'LineWidth',2)
    end
    plot([t(1) t(end)],[0 0],'Color',[.6 .6 .6]); 
    title(['elec' num2str(e) ': SSVEP CorrErr '])
    
end


% It is a good idea to check for differences in phase and amplitude together using a polar plot where we can plot the
% conditions as full complex numbers. Note that for electrode 9 you can see that there is a somewhat consistent phase
% lag for error relative to correct
for e=[Bestelectrode_SN Bestelectrode_CorrErr(4)] %[9 22]
    figure; polar(max(SV1corr(:))); hold on
    trange = find(t>1000 & t<1500);
    for c=1:4
        plot(mean(SV1corr(e,trange,c),2),'o','Color',colors(c,:),'LineWidth',2)
        plot(mean(SV1err(e,trange,c),2),'d','Color',colors(c,:),'LineWidth',2)
    end
    title(['elec ' num2str(e) ' Copmlex Corr/Err'])
    
end

