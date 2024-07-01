% Behavioural analysis of the continuous dots task
% We will just compare high to low coherence trials
clear
% where are the behavioural data files (saved from matlab script running the task):
datapath = 'H:\Rose\Behaviour\';

% define relevant trigger codes
targtrig = [101:104];
resptrig = [51 53];

% Target Triggers (coherent motion onset):
% 40% leftward target = 101;
% 40% rightward target =102;
% 80% leftward target = 103;
% 80% rightward target= 104;
% responses:
% 51 left button
% 53 right button

% Make handy vectors for transforming above trigger codes into individual codes for motion direction (1=left,
% 2=right) and coherence (1=40%, 2=80%):
trig2md = [1 2 1 2]; % we will subtract 100 from the triggers first
trig2coh = [1 1 2 2]; % similarly, translating trigger code to coherence

% list subjects:
allsubj = {'P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'RL01', };

% list different block types (even though we will ignore block type for purposes of this demo)
blockConds = {'HC' 'LC' 'LD' 'RD'}; % blocked conditions 
       
minRT = 0.2; % in sec; what is the minimum RT we would regard as a genuine reactino to the target? In this task with contunuous dots and slowly ramping targets 0.2 sec is a conservative minimum. Other tasks would be different - e.g. in fast discrete tasks we might accept anticipatory responses before evidence onset as legitimate responses - think about what's appropriate for YOUR task
maxRT = 2.5; % also in sec; maximum RT we consider legitimate (slow!) response to target, as opposed to a false alarm that occurred during ITI
% 2.5 sec is a reasonable maxRT here because althoghu targets were 1760 ms, subjects could be very delayed in
% responding. I chose this based on the RT histograms below.

histbin = [0:.1:2.5]; % bins for histograms of RTs
for s = 1:length(allsubj)
        
    % For each subject we're going to record the following for each trial each in its own vector:
    coh=[]; % coherence
    md=[]; % motion direction
    RT=[]; % RT
    rLR=[]; % response chosen - left=1, right=2
    
    for bc=1:length(blockConds)
        for b=1:4
            % Load mat file for current block:
            try, load([datapath allsubj{s} '/' allsubj{s} '_' blockConds{bc} num2str(b)]); catch; end
            % The data we need from these mat files are:
            % PTBevent: vector of trigger codes of the triggers that came up in each block
            % PTBeventT: vector of times (in sec since, I think, the beginning of the computer's life!)
            % corresponding to the same triggers
          
            % find the indices of the target triggers:
            targs = find(ismember(PTBevent,targtrig));

            % we're going to use local per-block versions of the 4 vectors:
            clear coh1 md1 RT1 rLR1
            % transform target triggers to respective coherence and motion codes:
            coh1 = trig2coh(PTBevent(targs)-100);
            md1 = trig2md(PTBevent(targs)-100);

            % get RTs
            for n=1:length(targs)
                RT1(n) = nan; % initialise. Let's code misses as nan
                rLR1(n) = nan;
                % find index of next response trigger:
                nextrespind = find(PTBeventT > PTBeventT(targs(n))+minRT & ismember(PTBevent,resptrig),1);
                thisRT = PTBeventT(nextrespind) - PTBeventT(targs(n)); % compute RT: response time minus target onset time
                if thisRT < maxRT
                    rLR1(n) = find(resptrig == PTBevent(nextrespind));
                    RT1(n) = thisRT;
                end
            end
            % append:
            coh = [coh coh1];
            md = [md md1];
            RT = [RT RT1];
            rLR = [rLR rLR1];
        end
    end

    % now get histograms for this subject for each coherence level:
    for c=1:2
        h(:,c,s) = hist(RT(find(coh==c & rLR==md)),histbin); % rLR==md means the chosen response matches the motion direction - i.e. we're just considering RTs of correct responses
    end
end
% plot RT histograms:
figure; hold on
plot(histbin,mean(h,3))
ylabel('number of trials')
ylabel('RT (sec)')