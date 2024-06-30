clear all
dbstop if error

% list the subjects:
allsubj = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10'};
blocks = {[1:10]
          [1:10]
          [1:10]
          [1:10]
          [1:10]
          [1:10]
          [1:10]
          [1:10]
          [1:10]
          [1:10]};

% Triggers
% fixation on - 4
% par.CD_GRATON = [1 2];
% par.CD_CUE = 6;
% par.CD_EVON = 100+[[1:numconds]; numconds+[1:numconds]];
% par.CD_EVOFF = 200;
% par.CD_GO = 3;
% par.CD_FB = 10;
% par.CD_BUTTONS = [12 13]; % Left and right mouse buttons

% all the relevant triggers that we'll be concerned with:
relevant_triggers = [4 1 2 6 101:110 200 3 10 12 13];

% define other relevant trigger codes
resptrig = [12 13];
stimtrig = [101:104];
anapar = struct();  % keep all analysis parameters in a structure that can be saved with epoch data

datafolder = 'C:\Users\Jacko\Desktop\College - 5th Year\Thesis\Jacks Analysis\Data\EMG\';  % Set path for Behaviour (.mat) data here
eegfolder = 'C:\Users\Jacko\Desktop\College - 5th Year\Thesis\Jacks Analysis\Data\EEG\';   % Set path for EEG (.bdf) data here
bigmatfolder = 'processingfiles\'; % Where the data will be stored

% EEG parameters:
fs = 512;   % sample rate
nchan=128;  % number of EEG channels
next=8;     % number of external channels (EOG etc)
% load in a structure 'chanlocs' containing the x,y,z locations of each of the 128 scalp channels in the cap.
load chanlocsBioSemi128; 
VEOGchans = [129 130]; % which channels are vertical EOG? 
HEOGchans = [131 132]; % horizontal at sides of eyes

anapar.LPF = 1;    % 1 = low-pass filter the data, 0=don't.
Chans2LPF = [1:nchan VEOGchans HEOGchans]; % what channels for lp?

% The low-pass filter, as previously created by Simon Kelly (2023)
fc=38;       % Low Pass Filter cutoff in Hz
wc = fc/fs*2*pi;
LPK=[];     % for low pass kernel
L=77;
for i=1:L
    n = i-ceil(L/2);
    if n==0, LPK(i)=wc/pi;
    else, LPK(i) = sin(wc*n)/(pi*n); end
end
LPK = LPK.*hamming(L)';
anapar.detrend_data = 1;     % detrend the data (1) or not (0)?

% define the epoch to extract:
epoch_limits_ms = [-700 3600]; 
ts = round(epoch_limits_ms(1)/1000*fs):round(epoch_limits_ms(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
t = ts*1000/fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds

%define the baseline window:
BLwin = [-188 0]; % multiple of SSVEP cycle

% Now loop through the subjects and extract the single trial ERPs and the other important informaiton about every trial:
for s = 1:length(allsubj)
    disp(['Extracting epochs for subject ' allsubj{s} '...'])

    % make empty vectors and matrices for each thing we want to save on every single trial:
    contrast = [];   % contrast used
    duration = [];   % duration of trial (deadline)
    condition = [];  % high or low contrast (in 1s and 2s format) 
    threshold = [];  % threshold required
    respLR = [];     % participant response, 1 = left, 2 = right
    corrLR = [];     % correct response, 1 = left, 2 = right
    rt = [];         % response time of trial (in seconds)
    lenBase = [];    % length of baseline data
    lenTrial = [];   % length of trial data
    lenPost = [];    % length of post EMG data
    totEMG = [];     % total EMG, concatenation of baseline, trial, post
    perf = [];       % performance of partcipant [in/correct, early, late]
    blocknum = [];   % keep track of block number so we can interpolate on 
                     % per-block basis (bad electrodes usually fixed between blocks)
    erp = [];        % the EEG epochs. Whereas all of the above are 1-D row vectors, this will have dimensions 136 channels x epoch timepoints x trial - so the 3rd dimension should match the length of all preceding vectors
    EEGRT = [];      % EEG rt; required further down
    EEGrespLR = [];  % EEG response; required further down
    firstTilt = [];  % for SSVEP
    
    % now list the eeg files to be used for each subject:
    filenames = [];  f=0; % initialize as empty a list of filenames and corresponding conditions
        for b=blocks{s}
            f=f+1;
            filenames{f} = [eegfolder '\'  allsubj{s} '\'  allsubj{s} '_B'  num2str(b) '.bdf'];
            matfilenames{f} = [datafolder '\' allsubj{s} '\'  allsubj{s} '_'  num2str(b) '.mat'];
        end

        for f=1:length(filenames)
        EEG = pop_biosig(filenames{f}); % read in EEG data, uses EEGLAB 
        numev = length(EEG.event); % total number of events
        clear trigs stimes
        
            for i=1:numev
                if length(EEG.event(i).type)>9
                trigs(i)= str2num(EEG.event(i).type(11:end));          
                else
                trigs(i)= str2num(EEG.event(i).type);   
                end
                stimes(i)=round(EEG.event(i).latency);  
            end
            
        % Shave off beginning and end bits based on triggers, because sometimes experimenter is delayed starting or
        % stopping the recording
        goodstuff_start = max(stimes(find(ismember(trigs,relevant_triggers),1)) - fs, 0); % go back 1 second before first trigger that matters
        goodstuff_end = min(stimes(find(ismember(trigs,relevant_triggers),1,'last')) + 2*fs, size(EEG.data,2)); % go 2 sec beyond the last relevant trigger
        EEG.data(:,goodstuff_end+1:end)=[];
        EEG.data(:,1:goodstuff_start)=[]; stimes = stimes-goodstuff_start;  % stimes is adjusted here also of course, so its sample points correspond to the now truncated continuous EEG
        
        % detrend?
        if anapar.detrend_data, EEG.data = detrend(EEG.data')'; end
        
        % LP Filter
        if anapar.LPF, for q=Chans2LPF, EEG.data(q,:)=conv(EEG.data(q,:),LPK,'same'); end; end % it's an FIR filter so we can implement directly by convolution
        
        % find indices of triggers corresponding to the events we want to time-lock to (targets)
        targs = find(trigs==6);

        [Contrast Duration Condition Threshold RespLR CorrLR Rt LenBase LenTrial LenPost TotEMG Perf FirstTilt] = get_trial_data(matfilenames{f});
        
        


    if length(Threshold)~=length(targs)
        disp('problem')
           dbstop
    end

      % Now loop through the single trials and grow the trial-descriptor vectors and 'erp' matrix
        for n=1:length(targs)
            if stimes(targs(n))+ts(end) > size(EEG.data,2), continue; end % first, check whether the epoch relative to the current event is even contained within the bounds
            blocknum = [blocknum; f];
            contrast = [contrast; Contrast(n)];   
            duration = [duration; Duration(n)];   
            condition = [condition; Condition(n)];  
            threshold = [threshold; Threshold(n)];  
            respLR = [respLR; RespLR(n)];     
            corrLR = [corrLR; CorrLR(n)];     
            rt = [rt; Rt(n)];         
            lenBase = [lenBase; LenBase(n)];    
            lenTrial = [lenTrial; LenTrial(n)];   
            lenPost = [lenPost; LenPost(n)];    
            totEMG = [totEMG; TotEMG(n)];     
            perf = [perf; Perf(n)];       
            firstTilt = [firstTilt; FirstTilt(n)]; % first tilt relevant
            
            % to SSVEP
            nextstimind = find(ismember(trigs,stimtrig) & stimes>stimes(targs(n)),1); % find the index of the next stimulus in trigs/stimes
            nextrespind = find(ismember(trigs,resptrig) & stimes>stimes(targs(n)),1); % find the index of the next response in trigs/stimes
            if ~isempty(nextrespind)
                EEGRT = [EEGRT stimes(nextrespind)-stimes(nextstimind)];
                EEGrespLR = [EEGrespLR trigs(nextrespind)-11]; % subtract 11 as a quick way to translate response trigger code to 1's and 2's for left/right
            else % if there WAS no next response, set the response parameters for this trial as 'not a number'
                EEGRT = [EEGRT nan];
                EEGrespLR = [EEGrespLR nan];
            end
            
            % Now extract the epoch
            ep = EEG.data(:,stimes(targs(n))+ts);
            %Baseline correction
            BLamp = mean(ep(:,find(t>BLwin(1) & t<=BLwin(2))),2);
            ep = ep - repmat(BLamp,[1,length(t)]);
            erp = cat(3,erp,ep); % now add the current epoch onto the growing 'erp' matrix by concatenation along the 3rd dimension  
        end
      
    end
    
    % Now that we have all trials from all blocks, save everything for this subject:
    save([bigmatfolder allsubj{s} '_raw'],'erp','t','contrast','duration','condition','threshold','respLR','corrLR','rt','totEMG','perf','blocknum','filenames','EEGRT', 'EEGrespLR', 'anapar');
end

