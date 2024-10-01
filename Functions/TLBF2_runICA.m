function [] = TLBF2_runICA(sub,exp,EEG)

%Lowpass filter eye channels; not done before - if noisy, ICA is not clean
% EEG1 = pop_eegfiltnew(EEG, 'locutoff', 1, 'plotfreqz',1,'channels',[]);

% EEG1 = EEG; 这里加1Hz滤波提高ICA质量
EEG1 = pop_basicfilter(EEG, 1:exp.nEEGchans, 'Cutoff', 1, 'Design', 'butter', 'Filter', 'highpass', 'Boundary', [], 'order', 4);

[EEG1] = pop_resample(EEG1,128); %Try downsampling just for ICA

EEG1 = EEG;
tic()
EEG1 = pop_runica(EEG1, 'icatype', 'runica', 'extended',0,'interrupt','on', 'chanind', exp.icaChans) % exp.icaChans);
toc()

EEG.icawinv = EEG1.icawinv;
EEG.icasphere = EEG1.icasphere;
EEG.icaweights = EEG1.icaweights;
EEG.icachansind = EEG1.icachansind;
% 
% data4pca = double(reshape(EEG1.data(1:exp.nEEGchans,:,:),[exp.nEEGchans,size(EEG1.data,2)*size(EEG1.data,3)]));
% %PCA 
% [pc,eigvec,sv] = runpca(data4pca, 10);

% Upload ICA components to full-sample dataset
EEG=  pop_saveset(EEG, ['ICA' EEG.filename], exp.filepath); %  for corrected


clear EEG;


end