function [EEG1,EEG2] = TLBF2_filter(sub,exp)


% Load preprocessed data and epoch around SL and responses
EEG = pop_loadset( [exp.filepath '\' exp.filterLab 'aac' exp.name num2str(sub) '.set'] );

EEG1 = pop_epoch( EEG, exp.trigg.SL, exp.stimEpoch, 'epochinfo', 'yes');
EEG1 = eeg_checkset( EEG1 );
filename=['SL_' EEG.filename];%EEG.filename = filename;
EEG1 = pop_saveset( EEG1,filename, exp.filepath);

% EEG2 = pop_epoch(EEG, {exp.trigg.response}, exp.respEpoch, 'epochinfo', 'yes');
% EEG2 = eeg_checkset( EEG2);
% filename=['RL_' EEG.filename];%EEG.filename = filename;
% EEG2 = pop_saveset( EEG2,filename, exp.filepath);

% EEG2 = pop_epoch(EEG, {'66'}, exp.respEpoch, 'epochinfo', 'yes');
% EEG2 = eeg_checkset( EEG2);
% filename=['EoL_' EEG.filename];%EEG.filename = filename;
% EEG2 = pop_saveset( EEG2,filename, exp.filepath);

EEG2 = pop_epoch(EEG, {'55'}, exp.respEpoch, 'epochinfo', 'yes');
EEG2 = eeg_checkset( EEG2);
filename=['RL_' EEG.filename];%EEG.filename = filename;
EEG2 = pop_saveset( EEG2,filename, exp.filepath);


% EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset
% 
% filename = ['e_' EEG.filename];
% EEG.filename = filename;
% 
% EEG = pop_saveset(EEG, filename, exp.filepath);
% % Display success message
% disp(['File  ' filename '  saved successfully.']);
end