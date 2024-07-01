% This script Goes through the big mats holding epoched data of all
% subjects and interpolates the bad channels listed in the bandchannel
% mats. This uses a function eeg_interp which is from EEGLAB, so make sure
% you have EEGLAB downloaded and included in the path.

clear all; close all

allsubj = {'P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16'}; 

bigmatfolder = 'bigmats/'; % where are the big mat files in which you saved the erps and trial info for each subject?
badchanfolder = 'badchannels/'; % and where is the 'ch2interp' for each subject? 
load chanlocsBioSemi128; % structure specifying the channel locations in 3d space - used for topography plotting etc
nchan = 128; % number of channels
next = 8;  % number of externals

for s = 9%1:length(allsubj)

    load([bigmatfolder allsubj{s} '_raw'])
    load([badchanfolder 'ch2interp' allsubj{s}])  % These came from running BadChannelCheck on the *raw.mats
              
    % EEG cleaning - interpolate bad channels identified in BadChannelCheck:
    for b=1:length(ch2interp) % for each block
        thisblock = find(blocknum==b); % identify trials of that block
        badchans=ch2interp{b};   % get the bad channels recorded during BadChannelCheck for that block
        if isempty(badchans), continue; end % don't go further if there were no bad channels
        load blankEEG;  % load a clean-slate EEG structure that has no data, but has all the right fields for EEGLAB to jazz with it
        EEG.nbchan=nchan; % set number of electrodes to just the number of scalp electrodes - we will not interpolate based on the externals!
        EEG.data = erp(1:nchan,:,thisblock); % feed it the epoched data for this block, again only the scalp channels in the cap
        EEG.pnts=length(t); % it seems to need this too - epoch length in sample pts
        EEG.trials=length(thisblock); % number of trials
        EEG.chanlocs = chanlocs(1:nchan); % it needs channel locations too
        EEG=eeg_interp(EEG,badchans,'spherical'); % this line does the actual interpolation
        erp(1:nchan,:,thisblock)=EEG.data; % now replace the relavant parts of the big 'erp' matrix with the interpolated version
        % Note the externals will still be sitting there in channels 129-136, unaltered.
    end
    % now re-save inerpolated version:
    save([bigmatfolder allsubj{s} '_intp'],'erp','t','modir','coh','cond','respLR','RT','blocknum','filenames','anapar');
end
