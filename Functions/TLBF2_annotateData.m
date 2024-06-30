function [] = TLBF2_annotateData(sub,exp)

% % Load converted data
EEG = pop_loadset([exp.filepath 'ac' exp.name num2str(sub) '.set']);



% Add channel locations
% Based on https://www.biosemi.com/headcap.htm. Manually edited for one
% participant, saved as biosemi_128_ext8.ced file.
load chanlocsBioSemi_128_ext8_noeyes_nosignal;

%empty string for 
elab = strings(exp.nEEGchans,1);

%write label for 1-128
for indChan = 1:exp.nEEGchans
    elab(indChan) = chanlocs(indChan).labels;
end
elab(end)='FCz'; clear ind*

EEG.chanlocs = chanlocs;
clear chanlocs;



EEG.data(137,:,:) = []; %Drop last channel;
EEG.nbchan = size(EEG.data,1)

% % Load metadata
% EEG = EP1_addMetadata(sub,exp, EEG, 0);
% 
% % Add trial idx to event file
% EEG = EP1_addTrialIdx(sub, exp, EEG,0);
% 
% % Check that triggers & metadata match 
% EEG = EP1_alignmentCheck(sub,exp,EEG);


EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset

filename = ['aac' exp.name num2str(sub)];
EEG = pop_saveset(EEG, filename, exp.filepath);
% Display success message
disp(['File  ' filename '  saved successfully.']);
end