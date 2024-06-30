function [] = TLBF2_convertData(sub,exp)

for block=1:8
    %Convert data from .bdf (BioSemi) format to .set (EEGlab) format
    EEG = pop_fileio([exp.database '/' exp.name num2str(sub) '_' num2str(block) '.bdf'], 'dataformat', 'auto');
    EEG = eeg_checkset( EEG );

    filename = ['c' exp.name num2str(sub) '_' num2str(block)];
    EEG = pop_saveset(EEG, filename, exp.database);
    % Display success message
    disp(['File  ' filename '  saved successfully.']);
end


end
