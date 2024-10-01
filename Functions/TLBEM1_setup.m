function [exp] = TLBF1_setup(exp); 

addpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\Functions');
cd('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2');

exp.name            = 'TL';

exp.behpath           = ['G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\EMG_bahivour\'];
exp.database          = ['G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Liang_EEGraw\']
exp.filepath          = ['G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\preprocessed_data\'];
exp.plotpath          = ['G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\plots\'];
exp.finalpath          = ['G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\finial_data\'];

exp.nsub              = length(exp.sub_id);
exp.fs                = 512 % AFTER downsampling

%% behaviour parameter
exp.blocknum=8;
exp.trialnum=128;

%% EEG parameter
% External channel info 
exp.chan.veog1 = 129; 
exp.chan.veog2 = 130; 
exp.chan.heog1 = 131; 
exp.chan.heog2 = 132; 
exp.chan.emgr1 = 133; 
exp.chan.emgr2 = 134; 
exp.chan.emgl1 = 135; 
exp.chan.emgl2 = 136; 

exp.nEEGchans = 128;
exp.nchans = 137; % Last channel is empty; 

exp.eegChans = [1:128];
exp.eogChans = [129,130,131,132];

exp.srate    = 512; %sampling rate 

% Processing parameters
% FILTERING 
exp.filter.lowerbound = 0.1; % When 0.1, files are labelled f01; when 0.01, files are just labelled f.
exp.filter.upperbound = 30;

if exp.filter.lowerbound == 0.01
    exp.filterLab = 'f';
elseif exp.filter.lowerbound == 0.1
    exp.filterLab = 'f01';
end

% EPOCHING
% Stim-locking params
exp.trigg.SL = {101, 102, 103, 104}; exp.stimEpoch = [-1.5, 3]; % 2s is the max. 
% Action-locking params 
exp.trigg.response = 10; exp.respEpoch= [-3.5, 1];
% Epoch labels
exp.epochs = {'SL', 'RL','EoL'};

% ICA channels
exp.icaChans = [1:128];
% Plotting parameters 
load chanlocsBioSemi_128_EOG4_Liang

exp.ImportantChans = [4 19 85 54 115]; % electrodes around CPz (where CPP is), Fz (CNV) and left and right motor cortex (see 'cap_128_layout_medium.jpg')
% Also mark front-most channels:
exp.FrontChans = [71 72 80 81 93 94 103 70 73 79 82 92 95 102]; % at the very front edge SD can be high not necessarily because channels are bad but because person is blinking /moving eyes. We don't want to interpolate GOOD electrodes that later will HELP us pick up artifacts that we want to throw out.
exp.rerefchan = [23]; % pick a channel to re-reference the data to and get a second picture of variance across channels - important just because sometimes the reference used to load the EEG might itself have been bad during the recording...

exp.cppChans = {'A1','A2', 'A3','A4', 'A5', 'A6','A18', 'A19', 'A20', 'A31', 'A32', 'B1', 'B2', 'B3', 'C1', 'D1', 'D15', 'D16'};
exp.cppCentral = {'A1','A2', 'A3', 'B1', 'B2','C1', 'D1', 'D15', 'D16'};
exp.cppParietal= {'A4', 'A5', 'A6','A18', 'A19', 'A20', 'A31', 'A32', 'B3'};

exp.cnvChans = {'C28', 'C18', 'C15', 'C27', 'C19', 'C14','C26','C20','C13'}; %'C25', 'C21', 'C12','C24', 'C22', 'C11'};
exp.MBLChans={'D11','D12','D13','D18','D19','D20','D17','D27','D28'};
exp.MBRChans={'B30','B31','B32','B21','B22','B23','B17','B18','B19'};


% ARTEFACT REJECTION 
exp.lowthresh = -150; 
exp.upthresh = 150;


for el = 1:length(exp.cppChans)
    thisCh = exp.cppChans(el);
    exp.cppIdx(el) = find(strcmp({chanlocs.labels}, thisCh));
end

for el = 1:length(exp.cppCentral)
    thisCh = exp.cppCentral(el);
    exp.cppCentralIdx(el) = find(strcmp({chanlocs.labels}, thisCh));
end

for el = 1:length(exp.cppParietal)
    thisCh = exp.cppParietal(el);
    exp.cppParietalIdx(el) = find(strcmp({chanlocs.labels}, thisCh));
end

for el = 1:length(exp.cnvChans)
    thisCh = exp.cnvChans(el);
    exp.cnvIdx(el) = find(strcmp({chanlocs.labels}, thisCh));
end

for el = 1:length(exp.MBLChans)
    thisCh = exp.MBLChans(el);
    exp.MBLIdx(el) = find(strcmp({chanlocs.labels}, thisCh));
end


for el = 1:length(exp.MBRChans)
    thisCh = exp.MBRChans(el);
    exp.MBRIdx(el) = find(strcmp({chanlocs.labels}, thisCh));
end