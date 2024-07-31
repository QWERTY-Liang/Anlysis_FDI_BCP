function [] = TLB1_AllBehaviour(exp)

% requiured output
% the behavioural data matrix should contain one row for each trial
% col 1: subject num
% col 2: block num
% col 3: contrast (high or low contrast)
% col 4: muscle (1=FDI, 2=BCP)
% col 5: trial outcome (1=correct, 2=error, 4=no response, 3 too early, 6 slow, 5 wrong muscle)
% col 6: RT in sec
% col 7: evshowtime =[]   % 0.2sec additional time after response
% col 8: participant response, 1 = left, 2 = right 0= on response
% col 9: correct response, 1 = left, 2 = right
% col 10: first tilt used for SSVEP
%%%%%% col xxx: total EMG trail not ready
%%%%%% col xxx: total EEG trail not ready


% make empty vectors and matrices for each thing we want to save on every single trial:
    subID=[];        % col 1: subject num
    blocknum = [];   % keep track of block number so we can interpolate on 
    contrast = [];   % high or low contrast  
    muscle_used= []; % col 4: muscle (1=FDI, 2=BCP)
    perf = [];       % performance of partcipant [in/correct, early, late 6 conditions in total
    % % % % % % % % % % par.CORRECT=1;
    % % % % % % % % % % par.INCORRECT=2;
    % % % % % % % % % % par.TOOEARLY=3;
    % % % % % % % % % % par.TOOLATE=4;
    % % % % % % % % % % par.Wrongmuscle=5;
    % % % % % % % % % % par.Slow=6;
    rt = [];         % response time of trial (in seconds)
    evshowtime =[]   % 0.2sec additional time after response
    totEMG = [];     % total EMG, concatenation of baseline, trial, post
     totEMG_bcp = [];
    respLR = [];     % participant response, 1 = left, 2 = right % used to check too slow
    corrLR = [];     % correct response, 1 = left, 2 = right % used to check too slow
    TiltSSVEP= [];   % first tilt used for SSVEP
    
    % duration = [];   % duration of trial (deadline)
    % threshold = [];  % threshold required
    % lenBase = [];    % length of baseline data
    % lenTrial = [];   % length of trial data
    % lenPost = [];    % length of post EMG data
    for sub = exp.sub_id(1:end)
        subID = [subID; sub*ones(exp.trialnum*exp.blocknum,1)];
        for block=1:exp.blocknum

            filename = [exp.name num2str(sub) '_' num2str(block) '.mat'];
            disp(['Extracting behaviour data for subject ' filename '...'])

            [Contrast Muscle Perf Rt EVtime TotEMG TotEMG_BCP RespLR CorrLR FirstTilt] = TLB1_get_trial_data(filename,exp);



            blocknum = [blocknum; block*ones(exp.trialnum,1)];
            contrast = [contrast; Contrast'];
            muscle_used= [muscle_used; Muscle'];
            perf = [perf; Perf'];
            rt = [rt; Rt'];         % response time of trial (in seconds)
            evshowtime =[evshowtime; EVtime'];   % 0.2sec additional time after response
            totEMG = [totEMG; TotEMG'];     % total EMG, concatenation of baseline, trial, post
            totEMG_bcp = [totEMG_bcp; TotEMG_BCP']
            respLR = [respLR; RespLR'];     % participant response, 1 = left, 2 = right % used to check too slow
            corrLR = [corrLR; CorrLR'];     % correct response, 1 = left, 2 = right % used to check too slow
            TiltSSVEP= [TiltSSVEP; FirstTilt'];   % first tilt used for SSVEP


        end
       
    end

% Create a matrix for the double variables
AllBehaviour= [subID, blocknum, contrast, muscle_used, perf, rt, evshowtime, respLR, corrLR, TiltSSVEP];

% Create variable names
    Varable_name={'subID' 'blocknum' 'contrast' 'muscle_used' 'perf' 'rt' 'evshowtime' 'respLR' 'corrLR' 'TiltSSVEP'};
    
    % Create tables with variable names
T_AllBehaviour = array2table(AllBehaviour, 'VariableNames', Varable_name);
    
    
 % save as .mat file   
    save([exp.behpath exp.name '_ALL_include_EMG' ],'AllBehaviour','totEMG','totEMG_bcp');
% Write to a CSV file
writetable(T_AllBehaviour, [exp.behpath exp.name '_ALL_ForR.csv']);
    disp(['Saved for' exp.name '_ALL_ for_R/ EMG' '.csv and .mat'])
end
