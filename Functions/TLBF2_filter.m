function [EEG] = TLBF2_filter(sub,exp)
%ver2: 27 Sep 2024
%1.加了50Hz notch filter
%2.低通降低以适应ICA


% % Load converted data
EEG = pop_loadset([exp.filepath 'aac' exp.name num2str(sub) '.set']);


% Filter
EEG  = pop_basicfilter( EEG,  1:exp.nEEGchans , 'Cutoff', [exp.filter.lowerbound  exp.filter.upperbound], 'Design', 'butter', 'Filter', 'bandpass', 'Boundary', [], 'order', 4);
% 50 Hz Notch Filter using simplenotch
EEG = pop_basicfilter(EEG, 1:136, 'Cutoff', [49 51], 'Design', 'butter', 'Filter', 'simplenotch', 'Boundary', [], 'order', 4);
EEG = eeg_checkset( EEG );
if exp.filter.lowerbound == 0.1
    filename = ['f01' EEG.filename];
elseif exp.filter.lowerbound == 0.01
    filename = ['f' EEG.filename];
end

EEG.filename = filename;



%EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset

%filename = ['aac' exp.name num2str(sub)];
EEG = pop_saveset(EEG, filename, exp.filepath);
% Display success message
disp(['File  ' filename '  saved successfully.']);
end