%%  5.1 Time frequency STFT
% Author: Liang Tong
% Date: 12/7/2024 

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
%% STFT prepare
fs=EEG.srate;%512
t=EEG.times;
erp=EEG.data;

epoch_limits_msTG = [-188*6 188*12];    % Target-locked epoch. We could use the same windows as were used for artifact rejection, or we might want different windows here - if so, don't forget that artifacts were only checked for in the window specified in IdentifyBadTrials!
tts = round(epoch_limits_msTG(1)/1000*fs):round(epoch_limits_msTG(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
tt = tts*1000/fs; % hence timebase in milliseconds, for plotting etc
% epoch_limits_msR =  [-188*3 188];    % Response-locked epoch; again using integer number of SSVEP cycles (1000/18.75 = 53.33)
% trs = round(epoch_limits_msR(1)/1000*fs):round(epoch_limits_msR(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
% tr = trs*1000/fs; % hence timebase in milliseconds



% Let's say we also want to compute spectral anmplitude as a function of time in the Mu/beta bands (reflects motor preparation when measured
% over motor cortices). One way to do such time-frequency analysis is the short-Time Fourier Transform - taking FFTs in a sliding window
% across time. First define parameters of this analysis:
fftlen = round(fs/21.5*15); % Window of how many sample points? If there is an SSVEP involved, whether or not you are interested in analyzing it, it is good to have all power related to the SSVEP isolated in a single frequency bin. This happens when you choose a window length that is an integer number of SSVEP cycles.
F = [0:fftlen-1]*fs/fftlen; % frequency scale, given window length (remember resolution = 1/window-duration)
ff = find((F>4 & F<21.45) | (F>21.55 & F<50)); % the indices of F that cover the spectral range/band of interest. Let's say we're interested in Mu and Beta bands combined. Note here I'm avoiding the SSVEP frequency (18.75hz in this example)
Ts = [-188*6:47:188*12]; % in msec, centered on what times do you want to measure spectral amplitude? i.e., where to center each consecutive window in time
% Tr = [-188*3:47:188]; % for response-locked

% Now in the rest of the code in the loop, we turn to computing and averaging Mu/Beta amplitude
% compute short time fourier transform (STFT) for each single trial:
STFT = []; % initialise
% for the stimulus-locked STFT, we can compute ffts across all trials at once - that's why there is no 'for n=1:ntr' loop.
for tt=1:length(Ts) % for each STFT timepoint (window centre) we'll compute the FFT for all trials at once
    [blah,samp] = min(abs(t-Ts(tt))); % find the sample point in the ERP epoch corresponding to the centre of the current FFT window
    spec = abs(fft(erp(:,samp-round(fftlen/2)+[1:fftlen],:),[],2))./(fftlen/2); % compute the magnitude of FFT (and scale it so that it's amplitude in same units as EEG
    % Save the result for each trial, just like the matrix 'erp'
    STFT(:,tt,:,:) = spec(:,ff,:);
end
% % response-locked (again processing each single trial):
% STFTr = nan(size(erp,1),length(Tr),ntr); % initalise. This time to NaN, because we want to retain same dimensions (a layer for each single trial) for indexing, but if R is such that we can't chop out a full STFT timecourse for a given trial, we mark that as NaN
% validrlockSTFT = zeros(1,ntr); % Again, make vector indicating whether valid for response locking (1) or not (0). The criteria for this are adjusted to account for windowing in STFT
% for n=1:ntr
%     [blah,RTsamp] = min(abs(t-RT(n)*1000)); % find the sample point closest to the RT for this trial. Remember RT is in sec
%     if RTsamp+Tr(1)*fs/1000-fftlen/2>0 & RTsamp+Tr(end)*fs/1000+fftlen/2<=length(t) & ~isnan(RT(n)) % Again, a trial is only validrlock=1 if the RT lies within the erp epoch you're using
%         validrlockSTFT(n)=1;
%         for tt=1:length(Tr)
%             spec = abs(fft(erp(:,RTsamp+round(Tr(tt)*fs/1000)-round(fftlen/2)+[1:fftlen],n),[],2))./(fftlen/2); % compute the magnitude of FFT (and scale it so that it's amplitude in same units as EEG
%             STFTr(:,tt,n) = mean(spec(:,ff),2);
%         end
%     end
% end
AllBehaviour_SL_bb=AllBehaviour_SL_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5.3 sorting left and right (FDI vs BCP)

selected_trials_1 = zeros(length(AllBehaviour_SL_bb),1);  % left C4-C3
selected_trials_2 = zeros(length(AllBehaviour_SL_bb),1);  % right C3-C4
count=0;
% switch%%%%%%%%%%calculate correct first
%     case 'CW_SL_b' %correct vs wrong; pre cue SL
        for i=1:length(AllBehaviour_SL_bb)
            %select trials
            if AllBehaviour_SL_bb(i,8)==1 && AllBehaviour_SL_bb(i,4)==1 &&AllBehaviour_SL_bb(i,9)==1%&& AllBehaviour_SL_bb(i,3)==0.07
                selected_trials_1(i)=1;
            elseif AllBehaviour_SL_bb(i,8)==2 && AllBehaviour_SL_bb(i,4)==1 &&AllBehaviour_SL_bb(i,9)==2%&& AllBehaviour_SL_bb(i,3)==0.07
                selected_trials_2(i)=1;
            end
            %exclude invalide(too early or wrong muscle)
            if AllBehaviour_SL_bb(i,5)==3 || AllBehaviour_SL_bb(i,5)==5
                selected_trials_1(i)=0;%change back to unselected
                selected_trials_2(i)=0;
                count=count+1;
            end
        end



        trl1=find(selected_trials_1==1);
        trl2=find(selected_trials_2==1);

%BCP
selected_trials_11 = zeros(length(AllBehaviour_SL_bb),1);  % left C4-C3
selected_trials_22 = zeros(length(AllBehaviour_SL_bb),1);  % right C3-C4
count=0;
% switch%%%%%%%%%%calculate correct first
%     case 'CW_SL_b' %correct vs wrong; pre cue SL
        for i=1:length(AllBehaviour_SL_bb)
            %select trials
            if AllBehaviour_SL_bb(i,8)==1 && AllBehaviour_SL_bb(i,4)==2 &&AllBehaviour_SL_bb(i,9)==1%&& AllBehaviour_SL_bb(i,3)==0.07
                selected_trials_11(i)=1;
            elseif AllBehaviour_SL_bb(i,8)==2 && AllBehaviour_SL_bb(i,4)==2 &&AllBehaviour_SL_bb(i,9)==2%&& AllBehaviour_SL_bb(i,3)==0.07
                selected_trials_22(i)=1;
            end
            %exclude invalide(too early or wrong muscle)
            if AllBehaviour_SL_bb(i,5)==3 || AllBehaviour_SL_bb(i,5)==5
                selected_trials_11(i)=0;%change back to unselected
                selected_trials_22(i)=0;
                count=count+1;
            end
        end



        trl11=find(selected_trials_11==1);
        trl22=find(selected_trials_22==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% compute
clear avMB1 avMB2 avMB11 avMB22
        avMB1(:,:,:) = mean(STFT(:,:,:,[trl1;trl2]),4);% correct-left
       % avMB2(:,:,:) = mean(STFT(:,:,:,trl2),4);% correct-right
        avMB11(:,:,:) = mean(STFT(:,:,:,[trl11;trl22] ),4);% wrong-left
        %avMB22(:,:,:) = mean(STFT(:,:,:,trl22),4);% wrong-right


      % ch = [115+1 54+1];
ch = [115 54];


      
       %% plot different condition same channel
       % Plot the spectrogram
       
figure (1);
imagesc(Ts, F(ff), ( squeeze(avMB1(ch(1),:,:)) -squeeze( avMB11(ch(1),:,:)) )' ); % Plot in decibels
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram - FDI-BCP left channel');
colorbar;
colormap jet;
% Set color limits
caxis([-0.5 0.5]); % Set color limits from -0.4 to 0.4
% Add vertical lines at specified time points
hold on;
xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s'});

% Add horizontal lines for mu and beta bands
yline(8, '--k', 'Mu lower');
yline(12, '--k', 'Mu-beta');
yline(30, '--k', 'Beta upper');

hold off;
% plot different channel

figure (2);
imagesc(Ts, F(ff), ( squeeze(avMB1(ch(2),:,:)) -squeeze( avMB11(ch(2),:,:)) )' ); % Plot in decibels
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram - FDI-BCP right channel ');
colorbar;
colormap jet;
% Set color limits
caxis([-0.5 0.5]);
% Add vertical lines at specified time points
hold on;
xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s'});
% Add horizontal lines for mu and beta bands
yline(8, '--k', 'Mu lower');
yline(12, '--k', 'Mu-beta');
yline(30, '--k', 'Beta upper');
hold off;

       %% plot different channel same condition
       % Plot the spectrogram
figure (1);
imagesc(Ts, F(ff), ( squeeze(avMB1(ch(1),:,:)) -squeeze( avMB1(ch(2),:,:)) )' ); % Plot in decibels
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram - FDI left-right');
colorbar;
colormap jet;
% Set color limits
caxis([-0.5 2]);
hold on;
xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s'});

% Add horizontal lines for mu and beta bands
yline(8, '--k', 'Mu lower');
yline(12, '--k', 'Mu-beta');
yline(30, '--k', 'Beta upper');

hold off;


% plot different channel

figure (2);
imagesc(Ts, F(ff), ( squeeze(avMB11(ch(1),:,:)) -squeeze( avMB11(ch(2),:,:)) )' ); % Plot in decibels
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram - BCP left-right');
colorbar;
colormap jet;
caxis([-0.5 2]);
hold on;
xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s'});
% Add horizontal lines for mu and beta bands
yline(8, '--k', 'Mu lower');
yline(12, '--k', 'Mu-beta');
yline(30, '--k', 'Beta upper');


hold off;

