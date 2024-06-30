function [] = TLBF2_appendData(sub,exp)


% Initialize an empty ALLEEG structure
ALLEEG = [];

for block=1:8

    % % Load converted data
    EEG = pop_loadset([exp.database 'c' exp.name num2str(sub) '_' num2str(block) '.set']);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, block); % Store the dataset in ALLEEG

end


% Merge datasets
EEG = pop_mergeset(ALLEEG, 1:block, 0);

EEG = eeg_checkset( EEG );% Check the consistency of the merged dataset

filename = ['ac' exp.name num2str(sub)];
EEG = pop_saveset(EEG, filename, exp.filepath);
% Display success message
disp(['File  ' filename '  saved successfully.']);

end
