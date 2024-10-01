function [EEG] = TLBF2_prereReference(sub,exp, EEG)

%EEG = pop_loadset([exp.filepath '/i' epoch '_' exp.filterLab 'aacTL' num2str(sub) '.set']);

% CAR reference (Common Average Reference - 128 electrodes)
EEG = pop_reref( EEG, [1:128] ,'exclude',[129:136] ,'keepref','on');
EEG = eeg_checkset( EEG );

% Bipolar re-reference
% See: https://eeglab.org/tutorials/ConceptsGuide/rereferencing_background.html

% Eye Channels (VEOG EXT 1-2); (HEOG EXT 3-4);
EEG.data(exp.chan.veog1,:) =  EEG.data(exp.chan.veog1,:) - EEG.data(exp.chan.veog2,:)
EEG.data(exp.chan.veog2,:) =  EEG.data(exp.chan.veog2,:) - EEG.data(exp.chan.veog1,:)
EEG.data(exp.chan.heog1,:) =  EEG.data(exp.chan.heog1,:) - EEG.data(exp.chan.heog2,:)
EEG.data(exp.chan.heog2,:) =  EEG.data(exp.chan.heog2,:) - EEG.data(exp.chan.heog1,:)

% % EMG channels (RIGHT HAND: EXT 5-6); (LEFT HAND: EXT 7-8);
% EEG.data(exp.chan.emgr1,:) =  EEG.data(exp.chan.emgr1,:) - EEG.data(exp.chan.emgr2,:)
% EEG.data(exp.chan.emgr2,:) =  EEG.data(exp.chan.emgr2,:) - EEG.data(exp.chan.emgr1,:)
% EEG.data(exp.chan.emgl1,:) =  EEG.data(exp.chan.emgl1,:) - EEG.data(exp.chan.emgl2,:)
% EEG.data(exp.chan.emgl2,:) =  EEG.data(exp.chan.emgl2,:) - EEG.data(exp.chan.emgl1,:)

% filename = ['r' EEG.filename];
% EEG.filename = filename;
   
end