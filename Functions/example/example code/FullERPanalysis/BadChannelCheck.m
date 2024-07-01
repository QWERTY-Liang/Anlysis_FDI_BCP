% This script opens matfiles with single-trial ERP epochs that have been extracted by ExtractEpoch.m,
% which were not yet re-referenced or had any channel interpolated ('raw') but may have been filtered/detrended, 
% and identifies "Bad" channels that should be marked for interpolation in the next step 
% NB This script does NOT do the trial-by-trial artifact rejection - that comes later. This step should only 
% flag channels that are bad due to poor connection/ faulty electrode - not ones that are 
% noisy due to blinks/eye movements (those are not bad electrodes - they reflect a bad subject, or possibly an artifact checking window that is too long - see AverageAndPlot stage)

% How do we tell when a channel is bad? There is not one correct way. This script implements Simon's way, which manually looks 
% at the variance of the signal on all channels and defines 'bad' as having abnormally high (or low!) variance.
% This procedure is going to seem quite subjective, but it makes you consider each subject/block individually and think through 
% any unexpected glitches you see, and how to salvage as many datasets as possible. 

% Note now that you can't reduce the risk to zero, of marking a channel here as 'bad' when it's not really, or failing to mark a truly
% bad channel as 'bad'. Rest assured that if the latter happens, you will see at the averaging step that lots of trial-rejection is happening 
% due to one or two offending electrodes, and you can always come back here and add those channels to the bad list and re-make the ch2interp for that subject
% and re-run the interpolation and CSD steps.

% An alternative approach is to use an Automated algorithm for bad channel detection like PREP, for which there is a function in 
% EEGLAB (google it to download this toolbox). You will need EEGLAB anyway for interpolation in the next step, and plotting topographies.
% My issue with automatic algorithms is that they are usually one-size-fits-all and tend to mark too little or too many channels as bad.
% But on the other hand they are far more efficient because you don't have to do so much eyeballing as here. However,
% whatever you decide to use, always LOOK AT THE DATA.
% For this manual approach, a critical principle is that you select channels for interpolation in an unbiased way - e.g.
% you don't interpolate more channels in one block type than another systematically. 
% Our experiments are typically repeated-measures, so it's not a big problem if one subject is noisier than another
% (which is the case anyway whether you use a blanket policy or not)

% cap_128_layout_medium.jpg will be useful here, to see where electrodes are on the head.

% Remember this example dataset is a particularly bad example, with lots of problems like bridged electrodes. The idea is that, having gone through
% this, you will be extra careful collecting YOUR data, and won't have as much heartache here at the analysis stage.

% We will run this script manually for each individual at a time. The top cell reads in the data and plots; then for each subject write notes 
% and a record of code used to explore bad channels in a separate cell below. This helps you keep a record to come back to if you need.
% So you'll run the first cell for subj = 'P01'; then make a 2nd cell to figure out and note the bad channels for P01, then come back and run cell
% 1 for P02, then run code in a 3rd cell for P02, then back to cell 1 for P03, and so on. It is all done for you in this example, but that isthe
% sequence in which I proceeded.
% for convenience list the subjects:
% 'P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'RL01', }; 
subj = 'P09';  % change this each time you want to run this cell on a new subject

bigmatfolder = 'bigmats/';
badchanfolder = 'badchannels/';

load chanlocsBioSemi128; % this loads a structure specifying the channel locations in 3d space - used for topography plotting etc
nchan = 128; % number of channels
next = 8;  % number of externals

% We're going to plot the standard deviation of each channel and spot ones that stick out. It's worthwhile to particularly
% scrutinize electrodes that might be particularly important in the final analysis.
ImportantChans = [4 19 85 54 115]; % electrodes around CPz (where CPP is), Fz (CNV) and left and right motor cortex (see 'cap_128_layout_medium.jpg')
% Also mark front-most channels:
FrontChans = [71 72 80 81 93 94 103]; % at the very front edge SD can be high not necessarily because channels are bad but because person is blinking /moving eyes. We don't want to interpolate GOOD electrodes that later will HELP us pick up artifacts that we want to throw out.
rerefchan = [23]; % pick a channel to re-reference the data to and get a second picture of variance across channels - important just because sometimes the reference used to load the EEG might itself have been bad during the recording...

% load in the epoched data
load([bigmatfolder subj '_raw']) 
    
% concatenate all epochs and get Standard deviation (SD) per channel
clear SD SDoz
for b=unique(blocknum) % for each block compute the SD
    trl=find(blocknum==b); % find the indices of the trials from the current block
    conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]); % concatenate all trials
    % Now one thing about these raw (as recorded) EEG signals is that the electrodes near the CMS-DRL reference electrodes will tend to be tiny just due to proximity to reference, and even when very noisy you'd miss them in this comparison of S.D.s for that reason.
    % So it's worth also looking at SD when referenced to somewhere else - electrode Oz (23) for example - just to make sure you haven't missed any bad channels
    conc2 = conc - repmat(conc(rerefchan,:),[nchan+next,1]); 
    for q=1:nchan+next % include externals here as well, because it might be a good idea to check whether those are noisy too
        SD(q,b) = std(conc(q,:)); % measure S.D. of each channel 
        SD2(q,b) = std(conc2(q,:));
    end
end
% are there any channels that stick out in terms of standard deviation?
% To check, plot the SD per channel, separately for each block (since often the electrodes are fixed with extra gel/poking around between blocks):
figure; 
subplot(2,1,1); hold on; plot(SD(1:nchan,:)); ylim([0 150]); ylabel('SD') % we only plot the channels in the cap because external electrodes are often higher variance (e.g. you might be recording EMG) and annoyingly set the scale so you always have to zoom in, and the purpose here is to identify channels for interpolation which is ALWAYS only the 128 cap channels
title(['subject ' subj])
% mark the important channels specified above - this is just a visual aid to know which you should particularly consider
for e=1:length(ImportantChans)
    plot([1 1]*ImportantChans(e),[0 max(SD(ImportantChans(e),:))],'k'); 
end
% mark the front edge channels as well, again a visual aid to know where blinks are likely to create higher variance (through no fault of the electrodes)
for e=1:length(FrontChans)
    plot([1 1]*FrontChans(e),[0 max(SD(FrontChans(e),:))],'r'); % front edge channels in RED
end
subplot(2,1,2); hold on; plot(SD2(1:nchan,:)); ylim([0 150]); ylabel('SD'); xlabel('electrode'); title(['Alternative reference elec ' num2str(rerefchan)])

numBL = length(unique(blocknum)) % number of blocks? (should be 16 each subject for this dataset)
ch2interp = cell(1,numBL); % empty cell array (repeated below but kept it here just in case forget in one of the cells below - each subject should start with blank slate)

% Another way to see if a channel that seems bad for certain blocks are actually bad channels, or is it a 
% bad subject doing something whacky on some trials, is to get SD per channel per trial:
SDct = squeeze(std(erp(1:nchan,:,:),[],2)); % after squeezing, it's the SD of each trial for 128 channels x numtrials
figure; imagesc(SDct,[0 50]);   % the last input of imagesc sets the color axis limits, and I chose them based on typical SD values from the first plot.
xlabel('trials (concatenated for all blocks)'); ylabel('electrode')
% Vertical stripes of high SD mean lots of channels are bad for a limited set of trials (so not a case for interpolation - just let artifact rejection criterion later exclude those trials) 
% Horizontal stripes of high SD means lots of trials having high SD for those specific channels - that's a
% case for interpolation

% General principles: In this step we're identifying bad CHANNELS for interpolation - channels that would be noisy even if the subject 
% were sitting perfectly still and not blinking, moving eyes, clenching, gritting, whatever. So beware of marking channels near the eyes 
% (e.g. 72, 80, 81, 93, 94) as bad when actually it's the SUBJECT who is bad in their blinking. We'll see examples of such subjects below.
% Often to make this distinction well, we have to additionally plot the raw time-course of individual channels.
% The other thing is that channels can be bad because they have suspiciously LOW variance too (SD close to zero) - this could happen because 
% it's not plugged in for example, and they should also be marked for interpolation.   

return; % this is here in case you mistakenly run the script rather than just the cell. 
% Below are cells where I write notes on each individual subject, and the code I use for extra plotting and channel-ordering, etc, you can copy &
% paste as you like.
%% Let's first take the example of P01:

ch2interp = cell(1,numBL); % makes an empty cell array for filling in the bad channels for this subject
% ch2interp will be saved for each subject and used in a later step to interpolate.

% This subject is one of the better ones.
% In the SD-vs-elec figure, some channels visibly stick out (use zoom in tool), e.g. 80. And this electrode stands out
% as a horizontal stripe in the surface plot. But this is one of the channels marked red, over the eyes, so it might be blinks rather than a bad channel.
% There is one block in SD-vs-elec figure that sticks out most - Let's find that block:
[sortedSD,I] = sort(SD(80,:),'descend') % This sorts the BLOCKS by descending variance for THAT channel. 
% Here I see that block 4 is the one with highest variance on channel 80,
% so let's view it by plotting the time-course of that electrode concatenated across trials:
b=4; % choose a block
trl=find(blocknum==b); % find the indices of the trials from that block
conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]); % concatenate all trials. conc has all channels.
% When checking a channel it's always good to plot it alongside others, to compare. In this case it seems from the SD-vs-electrode plot that for the same block 4 
% there are several channels with high SD that could be bad (all frontal), so plot them all, along with some neighboring ones that are not so high-variance (e.g. 70) to see what a 'better electrode' looks like:
elec2plot = [70 71 72 78 79 80 81 94 95 101 102 118 119]; 
figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot')); % offsets added to electrodes to tell them apart.
% Note that these are concatenated trials so there are dicontinuities (sudden jumps) at the transitions from one trial to the next
% It's not those jumps that are bad - they are not happening in the continuous EEG but rather are aritfacts of our concatenation - but when
% the jumps are big, it is because epochs are drifting far away from baseline, which is a sign of noisiness.

% Now which of these channels are "bad?"
% None of them are astoundingly bad compared to neighbours, and you can see that a large part of the high variance is just coming from blinks (zoom in to some of spikes)
% For the most part most of these channels look no worse than the low-variance electrode 70.
% You can see that some of the most severe drifting on certain trials happened to not one but a whole group of electrodes simultaneously 
% (and these elecs are all around the front (e.g. channel 80 = 3*32+16 is C16, right over the right eye))
% These are likely due to something like a movement (raised eyebrows? Eye-drooping?) and so the high variance in those cases are more to do with a bad artifact, 
% NOT bad electrodes, and it is better to let those trials be thrown out of the analysis at the later trial-averaging stage.
% You can see that electrode 80 shows some drifting of its own, a bit more than other electrodes, 
% so let's mark that for interpolation across all blocks:
for b=[4], % which block(s) to mark this channel for interpolation in? Just 4 in this case
    ch2interp{b}=[ch2interp{b} 80];  % adds '80' to the list of bad channels for block b
end

% Also, electrodes 118 and 119 are slightly more fuzzy than the others. This is quite typical of electrodes beside the ears. 
% Let's quickly check another block to see how general:
b=3; % the second noisiest for channel 80. Now all the same code as above:
trl=find(blocknum==b); % find trial indices
conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]); % concatenate
elec2plot = [70 71 72 78 79 80 81 94 95 101 102 118 119]; % choose elecs to plot
figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot')); % plot
% 80 doesn't stick out quite as much; 118 and 119 are again fuzzy but SD is not so big so let's leave them 

% After we're satisfied we've identified all bad channels (one elec in one block in this case), save them for this subject:
save([badchanfolder 'ch2interp' subj],'ch2interp')

%% Now subject P02 - go back and re-run the top cell with subj = 'P02'

ch2interp = cell(1,numBL); % empty cell array

% Also pretty ok
% Channel 17 looks like it sticks out for most blocks. Worst to best:
[sortedSD,I] = sort(SD(17,:),'descend')
% let's check it against neighbouring electrode that's equally far from CMS-DRL, say 21, for the worst block:
b=16; % choose a block
trl=find(blocknum==b); % find the indices of the trials from that block
conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]); % concatenate all trials
elec2plot = [17 21]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% that's definitely a bad electrode. List it for interpolation for every block:
for b=I(1:end), ch2interp{b}=[ch2interp{b} 17]; end 
% Look at the surface plot. Channel 17 does indeed stand out as a horizontal stripe of high variance. 
% You can see there is high variance in the 80s and 90s too, but not so localised to a single electrode, and of course we know that's where the blinks are picked up. 

% Electrode 22 looks bad for one block. 
% Let's plot it for its worst block along with the neighbour up from it (23). I've put all above code into one line for ease of running:
ch=22; [sortedSD,I] = sort(SD(ch,:),'descend'); b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]); elec2plot = [ch ch+1]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% looks like there was a 'pop' just on one trial, so I won't interpolate it.

% What about channel 40?
ch=40; [sortedSD,I] = sort(SD(ch,:),'descend'); b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]); elec2plot = [ch ch+1]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% Again seems to have just had that one pop. But in this case it seems to happen on MOST blocks (judging from surface plot + SD-vs-elec plots) 
% So let's also give it the chop:
for b=I(1:end), ch2interp{b}=[ch2interp{b} 40]; end 

% Check out channel 94:
ch=94; [sortedSD,I] = sort(SD(ch,:),'descend'); b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]); elec2plot = [ch ch+1]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% yet another isolated pop, and few blocks affected. Leave it. 

% And channel 102 seems to stick out in one block:
ch=102; [sortedSD,I] = sort(SD(ch,:),'descend'); b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]); elec2plot = [ch ch+1]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% Localised to one trial. Leave it.

% Channel 121?
ch=121; [sortedSD,I] = sort(SD(ch,:),'descend'); b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]); elec2plot = [ch ch+1]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% That looks like slightly more than an isolated pop. From the SD plot it looks like there are 4 blocks
% for which it sticks out, so mark for interpolation for those blocks:
for b=I(1:4), ch2interp{b}=[ch2interp{b} 121]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp') % now save

%% Now subject P03 - go back and re-run the top cell..

ch2interp = cell(1,numBL); % empty cell array

% This is poor data. Apart from the obvious noisy electrodes 34 and 69, there are many series of consecutive
% electrodes that seem to take almost the exact same SD value. This is because the signal is basically the
% same due to BRIDGING - the gel has made its way from one hole in the cap to at least one other, so they
% are electrically connected and form one circuit node. This happens most often when a subject has long
% straight hair - in those cases you must try to spread the hair as evenly as possible and use a tight
% enough cap that the electrode holders aren't floating far away from the scalp due to hair, and ensure
% that you squirt the gel carefully starting right at the scalp and retreating out as you gently eject the
% gel on the way out, forming a connection directly from scalp to electrode, but NOT with the next hole down
% This subject might not be useable as the activity of parietal and occipital areas of the scalp are
% effectively "short circuited." But you can keep going for now and make a note to check the signals and topographies
% of that subject on their own to check whether they are an outlier in terms of their signals. this must
% be done in an unbiased way - you shouldn't reject this subject just because their effects go in the wrong way. To avoid bias
% you can collapse across the conditions you hypothesise to be different, so you can't even see any such effect when assessing noisiness.

% Also note in this subject, electrodes 5 and 18 are practically silent. They are the two below CMS, and other electrodes 
% equally to CMS are NOT so low-variance, so I suspect there was also bridging from CMS to those two as well.
% I will mark some of the worst high-variance and flat electrodes for interpolation, but given the bridging, 
% this subject might need to be rejected entirely. It is worth proceeding and after epoch averaging, checking their signal topography 
% to see if there is anything sensible-looking about it, and decide then whether to reject. 
% largest variance elecs:
for b=1:numBL, ch2interp{b}=[ch2interp{b} 34]; end 
% for channel 69 it looks bad for just 4 blocks:
[sortedSD,I] = sort(SD(69,:),'descend');
for b=I(1:4), ch2interp{b}=[ch2interp{b} 69]; end
% and the flat electrodes:
for b=1:numBL, ch2interp{b}=[ch2interp{b} 5 18]; end 

% Let's look into electrodes 26 to 32:
[sortedSD,I] = sort(SD(29,:),'descend'); % picking 29 as the middle of the pack
% plot worst block, with some other lower-variance electrodes in the general vicinity:
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [9 26 29 32 39 40]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% There's a big artifact affecting all electrodes, but there is also much more drift in these bridged electrodes, so let's mark to interpolate:
for b=I(1:3), ch2interp{b}=[ch2interp{b} 26:32]; end % only 3 blocks really stick out, judging from SD-vs-elec plot.

% and 51 is high-var in one block:
[sortedSD,I] = sort(SD(51,:),'descend');
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [50 51 52]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% isolated incident -> let that trial be rejected for high amplitude later.

save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P04 

ch2interp = cell(1,numBL); % empty cell array

% has so many noisy channels that non-noisy ones are the minority. Surface plot shows basically all
% channels were noisy for most trials. Looking at concatenated epochs over any single block shows there's
% a saw-tooth shape for many electrodes, and this is because there is slow drift that is being chopped
% into individually-baseline-corrected epochs. Going back to the raw EEG file (ExtractEpochs.m), and looking
% at raw signals before detrending or filtering, you see some huge slow fluctuations happening at least
% once in most blocks. This could be the head slowly nodding as subject fell asleep, or perspiration, or
% something else. There is not much interpolation can do here, so I left this subject with no electrodes noted.
% They will probably have high artifact detection rate no matter what we do.
save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P05: 

ch2interp = cell(1,numBL); % empty cell array

% The first block is worthless - blue light was probably flashing and all

% Surface plot gives the best view: electrodes 32 and 34 both almost completely flat for most trials, and
% they are respectively below and above DRL, so I suspect the 3 were bridged. 35 (B3) was also noisy on most blocks.
for b=1:numBL, ch2interp{b}=[ch2interp{b} 32]; end 
for b=1:numBL, ch2interp{b}=[ch2interp{b} 34]; end 
% 35 also sticks out for most blocks
for b=1:numBL, ch2interp{b}=[ch2interp{b} 35]; end 
% more horizontal stripes in surface plot:
for b=1:numBL, ch2interp{b}=[ch2interp{b} 53]; end 
for b=1:numBL, ch2interp{b}=[ch2interp{b} 80]; end 
for b=1:numBL, ch2interp{b}=[ch2interp{b} 120]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P06

ch2interp = cell(1,numBL); % empty cell array
% The obvious 3:
for b=1:numBL, ch2interp{b}=[ch2interp{b} 30 31 32]; end 
% investigate spike at 81 - about 4 blocks seem to stick out, and bit of a stripe in surf plot:
[sortedSD,I] = sort(SD(81,:),'descend');
for b=I(1:4), ch2interp{b}=[ch2interp{b} 81]; end 
% what's going on at 94?
[sortedSD,I] = sort(SD(94,:),'descend');
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [92:96]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% definitely jumping around more than neighbours, so mark for the two worst blocks
for b=I(1:2), ch2interp{b}=[ch2interp{b} 94]; end 
% 95?
[sortedSD,I] = sort(SD(95,:),'descend');
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [92:96]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% again, noisier than neighbours, so:
for b=I(1:5), ch2interp{b}=[ch2interp{b} 95]; end % I'm choosing number of blocks to interpolate this channel based on how many SD spikes seem separated from the pack
% 99, 2 blocks seem high-variance compared to rest
[sortedSD,I] = sort(SD(99,:),'descend');
for b=I(1:2), ch2interp{b}=[ch2interp{b} 99]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P07

% Finally a subject with no really high SDs, though, zooming in, it looks
% like there was some bridging here.

ch2interp = cell(1,numBL); % empty cell array
save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P08

ch2interp = cell(1,numBL); % empty cell array

% 11 and 12 look high for two blocks
[sortedSD,I] = sort(SD(11,:),'descend');
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [11 12 13]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% isolated popping events that we probably want to let be rejected later

% But also notice that electrode 1 is zero-SD across the board (did you miss that? Don't forget to look for particularly LOW variance too!)
for b=1:numBL, ch2interp{b}=[ch2interp{b} 1]; end 
% .. and 5 and 18, which are just below CMS, look almost flat. 6, equally close, meanwhile, is much higher variance, so this must be
% bridging from CMS to 5 and 18.
for b=1:numBL, ch2interp{b}=[ch2interp{b} 5 18]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P09

ch2interp = cell(1,numBL); % empty cell array

% high SD around the front channels near eyes. Investigate, taking highest-SD block for elec 80 as example: 
[sortedSD,I] = sort(SD(80,:),'descend');
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [72 80 81 82 93]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% wicked blinker. Leave all those channels as they are clean (..ly showing the terrible artifacts of the subject!)

% investigate 101:
[sortedSD,I] = sort(SD(101,:),'descend');
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [100:105]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% a bit worse than neighbours. Let's just mark it for interp for that one block with highest SD:
for b=I(1), ch2interp{b}=[ch2interp{b} 101]; end 

% similarly:
[sortedSD,I] = sort(SD(104,:),'descend');
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [100:105]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% Doesn't stick out during the non-artifact periods. Leave it.

% I've commented out the next line because it's something I did AFTER identifyBadTrials - come back here when you come to it!
% for b=1:numBL, ch2interp{b}=[ch2interp{b} 69 93]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P10

% a lot of high-SD channels. Looking at surface plot, it seems there are
% just a certain set of trials coming up intermittently, where a whole
% bunch of channels are higher-SD, so this is probably a problem with the
% subject, not the electrodes. 46 and 58 seem especially bad so mark those
% also channel 1 is flat again.

ch2interp = cell(1,numBL); % empty cell array

for b=1:numBL, ch2interp{b}=[ch2interp{b} 1 46 58]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')


%% P11 

ch2interp = cell(1,numBL); % empty cell array

% In 3 of the blocks, SD was huge for all channels (except 1 which was
% flat), but the surface plot shows this comes from pretty isolated events.
% On opening the raw bdf files and plotting all channels, it was clear that there were some
% major but very short disaster events in the data, but often in between, the data were fine. 
% These are crazy artifact trials, not bad channels, so all that is
% required here is to interpolate the flat channel 1.

for b=1:numBL, ch2interp{b}=[ch2interp{b} 1]; end 

% 119 seems to stick out on some trials...
[sortedSD,I] = sort(SD(119,:),'descend');
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [119 120]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% no it looks the same as 120 when zoom in. That crazy stuff at the beginning must have caused that.

save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P12 

ch2interp = cell(1,numBL); % empty cell array

% nice signals but apparently a bad blinker.
% 113 sticks out for one block
[sortedSD,I] = sort(SD(113,:),'descend');
b=I(1); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [112 113 114]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% isolated event but lasts a while so:
for b=I(1), ch2interp{b}=[ch2interp{b} 113]; end 
% it's evident there and in SD-vs-elec plot that 112 is flat. 5, 17 and 18 also

for b=1:numBL, ch2interp{b}=[ch2interp{b} 5 17 18 112]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')


%% P13
ch2interp = cell(1,numBL); % empty cell array

% another blinker but no problematic-looking horizontal stripes. 1 is flat again:

for b=1:numBL, ch2interp{b}=[ch2interp{b} 1]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P14
ch2interp = cell(1,numBL); % empty cell array

% Bridging across many pairs of electrodes, but might be ok. Just flat at channel 1 again. It must have been plugged out because it was
% discovered to be faulty
for b=1:numBL, ch2interp{b}=[ch2interp{b} 1]; end 

% 103 sticks out for 3 blocks
[sortedSD,I] = sort(SD(103,:),'descend');
for b=I(1:3), ch2interp{b}=[ch2interp{b} 103]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P15
ch2interp = cell(1,numBL); % empty cell array

% 3 blocks with really terrible explosions (have to zoom right out in SD-vs-elec plot) but these are localised vertical stripes in surface
% plot. Let's look at horizontal stripes:
% 32:
[sortedSD,I] = sort(SD(32,:),'descend'); % not the top 3 we know are those explosions, so let's look at 4th worst:
b=I(4); trl=find(blocknum==b); conc = reshape(erp(:,:,trl),[nchan+next,length(trl)*length(t)]);
elec2plot = [31 32]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))
% That's pretty bad alright.
for b=I(4), ch2interp{b}=[ch2interp{b} 32]; end 

for b=1:numBL, ch2interp{b}=[ch2interp{b} 1]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')

%% P16 

ch2interp = cell(1,numBL); % empty cell array

% an example of a very good subject and recording!
% just the flat channel 1:
for b=1:numBL, ch2interp{b}=[ch2interp{b} 1]; end 

save([badchanfolder 'ch2interp' subj],'ch2interp')
%%

% I'll leave it at this 16 subjects and keep going with those, but as an exercise you should do the rest.

 % In general, if later when you compute and plot average ERPs you see very localised blotches in topographies, 
 % or you see lots of trials rejected due to one or two offending electrodes, you should come back here and check 
 % whether you might have missed that particular electrode, and add it to ch2interp and re-run Interp etc.

