% CONTINUOUS Dots task: This means dots are always on the screen in the block and every now and then become coherent
% dot motion stimuli are created in the style of Shadlen lab

% CONVENTION: Anything that is a general task parameter or a single parameter for the block is put in the "par" structure
% Everything that is actual DATA is saved as its own variable, e.g. RespT, etc

% NOTE: dots stimulus code is taken from Shadlen Dots and we therefore adopt parameters in vis deg * 10, as they do

try  % the whole task code is put in a 'try' so that if there is an error it can exit gracefully and tell you what's wrong, rather than crash the computer
    % ***************************************************** BASIC SET - UP
    clear all;
    close all; clc
        
    AssertOpenGL;  % Break and issue an error message if the installed Psychtoolbox is not based on OpenGL or Screen() is not working properly. From Shadlen dots code, not sure if necessary.
    commandwindow;  % puts focus on command window, so if any keyboard buttons are pressed, they'll come up in the comand window rather than in your program!
    
    % load the parameters that were saved on last run (if not there, then start up a "par" structure with defaults in a couple of fields)
    try, load parCDots, catch, par.subjID='XX'; par.blocktype='HC'; par.blocknum=0; par.recordEEG=0; par.useEL = 0; par.numtrials=32; par.practice=0; end % on first run of task, set defaults to put in dialog box
    
    % MONITOR STUFF:
    whichScreen = 0; % depending on operating system the monitor will be indexed differently - usually 0, 1 or 2 - have to try and see
    par.monitorwidth_cm = 32;   % monitor width in cm - measure your screen! 32cm in booth @ UCD 
    dist_cm = 57;  % viewing distance in cm
    par.videoFrate = 60; % Desired video frame rate (refresh rate) of the monitor - 75 Hz is typical
%     par.screenres = [1080 768]; % desired screen resolution [W, H]
   % par.screenres = [1920 1080]; % desired screen resolution [W, H]
    par.screenres = [1920 1080]; %[1024,768]
    Skipsynctests = 1;
    par.diodeflash = 0; % we sometimes flash a square in top left of screen for a photodiode to precisely track timing. No need in UCD.

    % Stimulus Parameters:
    par.fixationSize = 0.2; % the fixation point square size in DEGREES
    par.BGcolor = 0; % background color - black
    par.rgbFIX = [255 255 255];% white
    
    par.cohSet = [0.4 0.8]; % set of Coherences that will be randomly interleaved; fully coherent = 1
    par.dirSet = [180, 0];   %%% set of Dot directions (also randomly interleaved) - left, right
   
    
    par.correctresp = [1, 3]; % the correct response buttons corresponding to the dirSet
    par.respstrings = {'LEFT','MIDDLE','RIGHT'}; % these should correspond to the mouse you're using - some have a middle; if it doesn't delete the middle one and adjust correctresp accordingly
    
    %  TIMING (all in ms - in the task trial loop /1000 to sec) - if 75 Hz use multiples of 40 ms
    par.ITISet = [960 1600 2240 2880]; % the set of possible ITIs - inter-target intervals
    par.fixperiod = 600;  %% Time for just fixation cross at the beginning
    par.evdur =  1760; % ms target evidence duration in total (1760 ms = 33 SSMEP cycles), the first 480 ms (9 cycles); % evidence duration
    par.rampdur = 480; % ramp-up duration at start of target
    par.FPC = 4; % #Frames per SSMEP cycle!
    
    % PARAMETERS ABOUT RESPONSE TIME AND REWARDS: (a lot of them!)
    par.minRT = 0.2; % in sec - starting from the onset of the ramp
    par.maxRT = 1.9; % in sec - note it is good to extend a little beyond the offset of coherent motion to give them a chance to press slightly after it
    par.maxpts = 10; % cap the points per trial
    par.ptsDecRT = 0.5; % RT beyond which the reward drops linearly
    par.ptsDecSlope = 7; % how steeply does reward decline with RT
    par.MISSpenalty = -4; % miss alarm penalty
    par.FApenalty = -4; % false alarm penalty
    par.Errpenalty = -4; % error penalty

    % re-seed the random number generator so randomisation won't be same every time you start matlab
    par.rseed = sum(100*clock); % Seed random number generator
    rng(par.rseed,'v5uniform'); % v5 random generator
    
    % There are two parameter structures borrowed from Shadlen Dots code: screenInfo and dotInfo. Some settings:
    screenInfo.dontclear = 0; % 0 for clear drawing, 1 for incremental drawing (does not clear buffer after flip)
    dotInfo.numDotField = 1; % number of dot patches
    
    % all below settings are those used by Shadlen group typically
    dotInfo.apXYD = [0 0 50]; %  x, y coordinates, and diameter of aperture(s) in 10ths of visual degree (so 50 is 5 deg)
    dotInfo.speed = 50;  % vertical vectors, dots speed (10th deg/sec) for each dot patch
    % For the Eye_Tracker to check for fixation, we create an aperture and
    % define its X-Y coordinates:
    th = 0*pi/180;     % polar angle, from positive x axis (horizontal meridian) - in or converted to radians
    r = 8;              % in degrees - the eccentricity - absolute distance from the center of the screen. Decrease this if ppt keeps losing fixation
    [par.posx,par.posy] = pol2cart(th,r); % In degrees, x-y coordinates
    
    
    dotInfo.maxDotsPerFrame = 150;   % From Shadlendots.
    % The dots routine tries to maintain a constant dot density, regardless of aperture size.  However, it respects MaxDotsPerFrame as an upper bound.
    
    dotInfo.dotSize = 5; % dot size in pixels
    dotInfo.dotColor  = [255 255 255]; % white dots default
    
    if par.useEL % Eyelink
        el.FixWinSize = 1.0;
        par.checkELfixation = 1;
        eye_used=-1;
    end
    
    % EEG
    par.eeg_pPort = 'C010';              % Address for parallel port (To find this, go to device manager -> ports. Right click on the parallel port and select properties. Go to the "resources" tab and the address is the the first element of the "I/O range".
    par.eeg_sRate = 1024;                                                          % Sampling rate of EEG system in Hz
    par.trigdur = 2*ceil(1000/par.eeg_sRate)/1000; % to be sure triggers safely deliver, set pulse duration of at least 2 EEG samples
    
    %%PTB stuff
    par.verbosity = Screen('Preference','Verbosity',2);
    par.skipSyncTests = Screen('Preference', 'SkipSyncTests',  Skipsynctests);
    par.visdeb = Screen('Preference', 'VisualdebugLevel', 3);
    
    [par.scresw, par.scresh]=Screen('WindowSize',whichScreen);  % Get screen resolution
    center = [par.scresw par.scresh]/2;     % useful to have the pixel coordinates of the very center of the screen (usually where you have someone fixate)
    par.hz=Screen('FrameRate', whichScreen,1);
    disp(['SCREEN RESOLUTION ' num2str(par.scresw) ' x ' num2str(par.scresh)])
    disp(['MONITOR REFRESH RATE ' num2str(par.hz) ' Hz'])
    
    if ~(par.scresw==par.screenres(1) & par.scresh==par.screenres(2)); disp(['Wrong resolution - set monitor to ' num2str(par.screenres(1)) ' x ' num2str(par.screenres(2))]); return; end
    
    if abs(par.hz-par.videoFrate)>1
        error(['The monitor is NOT SET to the desired frame rate of ' num2str(par.videoFrate) ' Hz. Change it.'])
    end
    
    % Get devices
    dev.mice = sort(GetMouseIndices);
    dev.mouse = dev.mice(1); % this should be the primary mouse - need to check though
    dev.mouseExperimenter = max(dev.mice);
    
    PsychHID('KbQueueCreate',dev.mouse);         % Create an event queue for the mouse (to listen out for mouse events)
    PsychHID('KbQueueStart',dev.mouse);            % Start listening out for mouse events
        
    cm2px = par.scresw/par.monitorwidth_cm;  % multiplication factor to convert cm to pixels
    deg2px = dist_cm*cm2px*pi/180;      % multiplication factor to convert degrees to pixels (uses aproximation tanT ~= T).
    
    % Rectangles are specified by [lower-left-xcoord lower-left-ycoord upper-right-xcoord upper-right-ycoord]. For fixation dot:
    fixSizePx = round(par.fixationSize*deg2px);
    fixRect = [center-floor(fixSizePx/2) center+(fixSizePx-floor(fixSizePx/2)-1)];
    
    dlg_title = 'Dot Direction Discrimination Task';
    while 1
        prompt = {'Enter Subject ID','Block type (which coh.dir more prob? LC, HC, LD, RD)','blocknumber (default below is last saved)','EEG? (1=yes, 0=no)','Eye tracker? (1=yes, 0=no)','Number of trials (multiple of numcoh x numdir)', 'Practice? (1=yes, 0=no)'};
        def = {par.subjID, par.blocktype, num2str(par.blocknum), num2str(par.recordEEG), num2str(par.useEL), num2str(par.numtrials),num2str(par.practice)};
        answer = inputdlg(prompt,dlg_title,1,def);
        par.subjID = answer{1};
        par.blocktype = answer{2};
        par.blocknum = str2num(answer{3});
        par.runID = [answer{1} '_' answer{2} answer{3}];
        par.recordEEG = str2num(answer{4});
        par.useEL = str2num(answer{5});
        par.numtrials = str2num(answer{6});
        par.practice = str2num(answer{7});
        
        if exist([par.runID '.mat'],'file') 
            dlg_title = [par.runID '.mat EXISTS ALREADY - CHOOSE ANOTHER, OR DELETE THAT ONE IF IT IS RUBBISH'];
        else
            break;
        end
    end
    
    % *********************************** BASIC SET - UP
    
%%%%%%%%%%%%%%% In case is a practice block, lower the coherence gradually:
%     Also make a prediction of their performance in Low Coh blocks. This
%     is done further on the script
 %%%%%%%%%%%%%%%%%
    if par.practice
            prompt2 = {'Enter coherences for the block: Hard and Easy'}; 
            def2 = {num2str(par.cohSet)};
            answer2 = inputdlg(prompt2,dlg_title,1,def2);
            par.cohSet = str2num(answer2{1}); % set of Coherences that will be randomly interleaved; fully coherent = 1
    end
 %%%%%%%%%%%%%%%%%        
    if par.useEL,
         par.edfFile=[par.runID(1:min(length(par.runID),8)) '.edf']; %make sure naming isnt too long that it cant be saved by EL
%         % calibrate:
        ELCalibrateDialog
    end
    %ListenChar(2); % Disable key output to Matlab window
    if par.recordEEG
        % Wake up the port
        eeg.obj = io64;                                   % obj = io32;
        eeg.status = io64(eeg.obj);               % stat = io32(obj);
        
        % Address to the port the way it likes it..
        port = hex2dec(par.eeg_pPort);
    end
   
    
%%%%%%%%%%%%%%% In case I want to debug while the program keeps running (openScreen):
%     PsychDebugWindowConfiguration
 %%%%%%%%%%%%%%%%% 
    
    % Open a graphics window on the main monitor
    [screenInfo.curWindow, screenInfo.screenRect] = Screen('OpenWindow', whichScreen, par.BGcolor);
    window = screenInfo.curWindow; % just to use a shorter word!
    Screen('TextSize', window, 18);
    
    % Get the refresh rate of the screen (seems redundant but using Shadlen code here - maybe this is more precise)
    screenInfo.frameDur = Screen('GetFlipInterval',screenInfo.curWindow);
    screenInfo.monRefresh = 1 / screenInfo.frameDur;
    
    % Get screen center coordinates (pixels)
    screenInfo.center = [screenInfo.screenRect(3), screenInfo.screenRect(4)]/2;   %%same as center above
    
    % Determine pixels per degree (ppd), which unit is derived as
    % (pix/screen) * (screen/rad) * (rad/deg) = (pixels per degree)
    screenInfo.ppd = pi * screenInfo.screenRect(3) / atan(par.monitorwidth_cm/dist_cm/2) / 360;
    % and more dots parameters...:
    apD = dotInfo.apXYD(:,3); % diameter of aperture - note still in 10ths of deg
    % Note the factor of 10 below comes in because Shadlen lab uses Rex which needs integers
    % find dot field centers in pixels, start by placing at the screen center... (Note in PTB y is inverted - pos on bottom, neg. on top)
    dotfieldcenter = repmat(screenInfo.center,dotInfo.numDotField); % a row for every dot field
    dotfieldcenter = [dotfieldcenter(:,1) + dotInfo.apXYD(:,1)/10*screenInfo.ppd, ... 
                      dotfieldcenter(:,2) - dotInfo.apXYD(:,2)/10*screenInfo.ppd]; % where you want the center of the aperture
    dotfieldcenter(:,3) = dotInfo.apXYD(:,3)/2/10*screenInfo.ppd; % add diameter (SK: RADIUS??)
    % Note this is called "center" in Shadlen code but we use that term for screen center so we're renaming it dotfield center here
    d_ppd = floor(apD/10 * screenInfo.ppd);	% diameter of aperture in pixels

    % ndots is the number of dots shown per video frame. Dots will be placed in a square of the size of aperture.
    % - Size of aperture = Apd*Apd/100  sq deg
    % - Number of dots per video frame = 16.7 dots per sq deg/sec,
    % When rounding up, do not exceed the number of dots that can be plotted in a video frame (dotInfo.maxDotsPerFrame).         
    ndots = min(dotInfo.maxDotsPerFrame, ...
        ceil(16.7 * apD .* apD * 0.01 / screenInfo.monRefresh));  % SK: what is the 16.7 I wonder? Should our number be different if we have 75Hz refresh rate? 
    
    %  ************************************************* CODES AND TRIAL SEQUENCE
    % Triggers are the numbers that will appear in the EEG on certain events. They can be whatever you want, as long as you
    % know them when you analyze the EEG. We have a convention of using 4 for fixation on...
    par.CD_StartTask = 1;
    par.CD_FIXON = 4;
    par.CD_DOTSON = 5;
    par.CD_TG = reshape(100+[1:length(par.dirSet)*length(par.cohSet)*dotInfo.numDotField],[length(par.dirSet),length(par.cohSet),dotInfo.numDotField]); % 
    par.CD_BUTTONS = [51 53]; % [left-button right-button]. 
    
    numcohs = length(par.cohSet);
    
    switch par.blocktype
        case 'HC' % high coherence most often
            par.DotDirProb= [50 50];
            par.DotCohProb = [25 75];
        case 'LC' % low coherence most often
            par.DotDirProb= [50 50]; 
            par.DotCohProb = [75 25];
        case 'LD' % Left direction most often
            par.DotDirProb= [75 25]; 
            par.DotCohProb = [50 50];
        case 'RD' % Right direction most often
            par.DotDirProb= [25 75]; 
            par.DotCohProb = [50 50];
    end
    
    if par.DotDirProb(1)~=par.DotDirProb(2) % biased direction, then...
        tgcoh=repmat([1:2],1,par.numtrials/2); % coherence, balanced. 1 is low coh and 2 is high coh
        tgdir=[ones(1,par.numtrials*par.DotDirProb(1)/100) 2*ones(1,par.numtrials*par.DotDirProb(2)/100)];
    elseif par.DotCohProb(1)~=par.DotCohProb(2) % biased coherence, then...
        tgdir=repmat([1:2],1,par.numtrials/2); % coherence, balanced
        tgcoh=[ones(1,par.numtrials*par.DotCohProb(1)/100) 2*ones(1,par.numtrials*par.DotCohProb(2)/100)];
    end
%     
%     % Make a sequence of trial conditions where every crossed condition is represented in the right proportion in the block:
%     tgdir = repmat([1:length(par.dirSet)],[1,floor(par.numtrials/length(par.dirSet))]); 
%     tgcoh=[]; for c=1:length(par.cohSet), tgcoh = [tgcoh c*ones(1,length(par.dirSet))]; end
%     tgcoh = repmat(tgcoh,[1,floor(par.numtrials/(length(par.dirSet)*length(par.cohSet)))]);  % index of coherence level on each trial 1,2, etc
%     % tgside only applies when there is more than one dot patch - on left and right for example
    tgside=[]; for d=1:dotInfo.numDotField, tgside = [tgside d*ones(1,length(par.dirSet)*length(par.cohSet))]; end
    tgside = repmat(tgside,[1,floor(par.numtrials/(length(par.dirSet)*length(par.cohSet)*dotInfo.numDotField))]); % 1 = targ on left, 2 = targ on right
        
    % randomize the ordering:
    r = randperm(par.numtrials);
    tgdir = tgdir(r);
    tgcoh = tgcoh(r);
    tgside = tgside(r);  % leaving this in in case we ever want two patches in future
    
    if par.useEL
        %%%%%%%%% EYETRACKING PARAMETERS
        par.FixWinSize = 4;    % RADIUS of fixation (circular) window in degrees
        par.TgWinSize = 3;    % RADIUS of fixation (circular) window in degrees
        ELsetupCalib
        Eyelink('Command', 'clear_screen 0')
        Eyelink('Command', 'draw_box %d %d %d %d 15', floor(center(1)-deg2px*par.FixWinSize), floor(center(2)-deg2px*par.FixWinSize), floor(center(1)+deg2px*par.FixWinSize), floor(center(2)+deg2px*par.FixWinSize));
        Eyelink('Command', 'draw_box %d %d %d %d 15', floor(center(1)+deg2px*par.posx-deg2px*par.TgWinSize), floor(center(2)-deg2px*par.posy-deg2px*par.TgWinSize), floor(center(1)+deg2px*par.posx+deg2px*par.TgWinSize), floor(center(2)-deg2px*par.posy+deg2px*par.TgWinSize));
        Eyelink('Command', 'draw_box %d %d %d %d 15', floor(center(1)+deg2px*-par.posx-deg2px*par.TgWinSize), floor(center(2)-deg2px*par.posy-deg2px*par.TgWinSize), floor(center(1)+deg2px*-par.posx+deg2px*par.TgWinSize), floor(center(2)-deg2px*par.posy+deg2px*par.TgWinSize));
        % figure out which eye is being tracked:
        if eye_used == -1
            eye_used = Eyelink('EyeAvailable');
        end
    end

    % *************************************************** START TASK
    % Present Instructions:
    HideCursor
    txtstart = 0.04*par.scresw;
    Screen('TextSize', window, 20)
    if dotInfo.numDotField==1
        Screen('DrawText', window, ['For the whole block, you will see a field of dots'], txtstart, 0.1*par.scresh, 255);
    elseif dotInfo.numDotField==2
        Screen('DrawText', window, ['For the whole block, you will see two fields of dots'], txtstart, 0.1*par.scresh, 255);
    end
    Screen('DrawText', window, 'At unpredictable times (roughly every 1-4 sec), the dots ', txtstart, 0.15*par.scresh, 255);
    Screen('DrawText', window, 'will move coherently to the left or right.', txtstart, 0.2*par.scresh, 255);
    Screen('DrawText', window, 'You must decide the direction of the dots each time, and', txtstart, 0.25*par.scresh, 255);
    Screen('DrawText', window, ['Press the ' par.respstrings{par.correctresp(1)} ' mouse button with your ' par.respstrings{par.correctresp(1)} ' hand if motion'], txtstart, 0.3*par.scresh, 255);
    Screen('DrawText', window, 'is LEFT, or', txtstart, 0.35*par.scresh, 255);
    Screen('DrawText', window, ['Press the ' par.respstrings{par.correctresp(2)} ' mouse button with your ' par.respstrings{par.correctresp(2)} ' hand if motion'], txtstart, 0.4*par.scresh, 255);
    Screen('DrawText', window, 'is RIGHT', txtstart, 0.45*par.scresh, 255);
    Screen('DrawText', window, 'The faster you respond, the more points you earn for each', txtstart, 0.5*par.scresh, 255);
    Screen('DrawText', window, 'correct target.', txtstart, 0.55*par.scresh, 255);
   % Screen('DrawText', window, 'You will lose points when you miss a target, click the wrong button, or when you click when there is no target', txtstart, 0.55*par.scresh, 255);
    Screen('DrawText', window, 'Stay fixated on central dot for the duration of the block, ', txtstart, 0.6*par.scresh, 255);
    Screen('DrawText', window, 'and if you need to blink, try to time it to just after you ', txtstart, 0.65*par.scresh, 255);
    Screen('DrawText', window, 'have clicked to a target', txtstart, 0.7*par.scresh, 255);
    
    switch par.blocktype
       case 'HC' % high coherence most often
             DrawFormattedText(window,'This is an EASY block.', txtstart, 0.8*par.scresh, [0, 255, 255])
%             Screen('DrawText', window, 'This is an EASY block.', txtstart, 0.8*par.scresh, 134);
       case 'LC' % low coherence most often
            Screen('DrawText', window, 'This is a HARD block.', txtstart, 0.8*par.scresh, [0, 255, 255]);
       case 'LD' % left direction most often
            Screen('DrawText', window, 'In this block LEFT motion occurs for every 3/4 trials.', txtstart, 0.8*par.scresh, [0, 255, 255]);
       case 'RD' % left direction most often
            Screen('DrawText', window, 'In this block RIGHT motion occurs for every 3/4 trials.', txtstart, 0.8*par.scresh,[0, 255, 255]);
    end
    
    
    Screen('DrawText', window, 'Press to begin', txtstart, 0.90*par.scresh, 255);
    Screen('Flip', window);     % This Flip function is important - you can draw all kinds of things on the "window" but they won't appear to the subject on the screen until the Flip.
     
    clear FixonT dotsonT cohstartT RespLR RT DotDir ButtT
    
    % Initialize some things that we'll record and save on a continual or
    % trial by trial basis (clear any lingering incarnations):
    PTBevent = [];  % PsychToolBox event code
    PTBeventT = []; % time at which each event occurred
    nPTBevent= 0; % number of events
    FixonT= nan(1,par.numtrials); %Time fixation
    dotsonT= nan(1,par.numtrials); %Dots on screen time
    cohstartT = nan(1,par.numtrials); %Target onscreen time
    RespLR = nan(1,par.numtrials); % on every trial, what's the response choice (first click of the trial), 1=left, 3=right
    RT = nan(1,par.numtrials); % on each trial, RT in sec
    Dotdir= nan(1,par.numtrials); % direction of the coherent motion/correct answer
    ButtT= nan(1,par.numtrials); %Button press time for  Response trigger
    
    % Wait for the user to press a button to begin task:
    [clicks,x,y,whichButton] = GetClicks(whichScreen,0); % GetClicks is a PTB function that 
                               % waits for an action on the mouse.
    StartTime = GetSecs;
                                    
    PsychHID('KbQueueFlush' ,dev.mouse,2)
    RespT(1) = GetSecs;
    Resp(1)=whichButton;    % The first response will be the one that sets the task going, after subject reads instructions
    
    %%%%%%%%%%%%%%%%%%%% START TRIALS
        
    % initial lead-in:
    Screen('FillRect',window, par.BGcolor); % screen blank
    Screen('Flip', window);
    WaitSecs(1);
    if par.recordEEG
        io64(eeg.obj,port,par.CD_StartTask);
        WaitSecs(par.trigdur);
        io64(eeg.obj,port,0);
    end
    if par.useEL
        Eyelink('Message', 'TASK_START');
    end
    
    saved_dots = cell(par.numtrials,dotInfo.numDotField); % initialise saved dots
    
    % Initial Fixation interval before dots start:
    Screen('FillRect',window, par.rgbFIX, fixRect);

    [VBLTimestamp, FixonT] = Screen('Flip', window, [], 0);  % "Dont clear" = 1 means it won't clear the framebuffer after flip (will  keep fixation point)
    if par.recordEEG,  io64(eeg.obj,port,par.CD_FIXON); WaitSecs(par.trigdur);  io64(eeg.obj,port,0); end
  % If using eyelink, Wait until Gaze enters fixation window to start trial... (this was adapted from EyelinkExample.m)
    if par.useEL
        Eyelink('Message', ['FIX_ON' num2str(par.CD_FIXON)]);
        while 1
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                evt = Eyelink( 'NewestFloatSample');
                % get current gaze position from sample
                x = evt.gx(eye_used+1)-center(1); % +1 as we're accessing MATLAB array
                y = -(evt.gy(eye_used+1)-center(2));
                % do we have valid data and is the pupil visible?
                if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                    if sqrt(x^2+y^2)<deg2px*par.FixWinSize
                        break;
                    end
                end
            end
        end
    end
    
    nPTBevent = nPTBevent+1;    % increment the number of events
    PTBevent(nPTBevent) = par.CD_FIXON;    % record the event that just happened (using the arbitrarily defined codes above)
    PTBeventT(nPTBevent) = FixonT;
   
    WaitSecs(par.fixperiod/1000 - .9*1/par.hz); % most of a refresh chopped off wait time to make sure it catches refresh it's supposed to
   
    for n=1:par.numtrials
        disp(['TRIAL ' num2str(n)])
        quit=0;

        % Set up trial
        dotInfo.dir = par.dirSet(tgdir(n));
        iti = par.ITISet(ceil(rand*length(par.ITISet))) -1; % -1 ms is there because otherwise it will probably have to wait until the NEXT refresh, so it will be out by one whole refresh...

        numFrames = round((par.evdur + iti)/1000*par.videoFrate); % total number of dots frame in this trial (ITI plus target)
        EvidenceOnsetFrame = round(iti/1000*par.videoFrate)+1;
        numRampFr = round(par.rampdur/1000*par.videoFrate);

        FixBreak=0; RespMade=0;  StimPresented=-1; dotsonT(n)=0; 

        displayed_dots = cell(numFrames,dotInfo.numDotField);

        %%Set up dots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Don't worry about pre-allocating, the number of dot fields should never be large enough to cause memory problems.
        % SK: the following seems to come up with random positions within square of [0 0 1 1] ("ss") for 3 sets of dots and
        % jump size assuming they're all coherent. These are dummy "last positions" that the dynamic code below can work
        % from, i.e. these ones are never actually plotted, they're just for initialisation.
        for df = 1 : dotInfo.numDotField
            % dxdy is an ndots x 2 matrix that gives jumpsize in units on 0..1
            %   deg/sec * ap-unit/deg * sec/jump = ap-unit/jump
            dxdy{df} = repmat((dotInfo.speed(df)/10) * (10/apD(df)) * ...
                (3/screenInfo.monRefresh) * [cos(pi*dotInfo.dir(df)/180.0), ...
                -sin(pi*dotInfo.dir(df)/180.0)], ndots(df),1);
            ss{df} = rand(ndots(df)*3, 2); % array of dot positions raw [x,y]
            % Divide dots into three sets
            Ls{df} = cumsum(ones(ndots(df),3)) + repmat([0 ndots(df) ndots(df)*2], ndots(df), 1);
            loopi(df) = 1; % loops through the indices of the three sets of dots - 1,2,3,1,2,3,1,2,3,1,...
        end

        dontclear = screenInfo.dontclear;

        for frame = 0:numFrames % play the whole duration no matter the RT 
            if 0%par.useEL  % want to check whether fixation is broken and quit trial if so?
%                if par.AbortOnFixationBreak
                    if Eyelink( 'NewFloatSampleAvailable') > 0
                        % get the sample in the form of an event structure
                        evt = Eyelink( 'NewestFloatSample');
                        % get current gaze position from sample
                        x = evt.gx(eye_used+1)-center(1); % +1 as we're accessing MATLAB array
                        y = -(evt.gy(eye_used+1)-center(2));
                        if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                            if sqrt(x^2+y^2)>deg2px*par.FixWinSize
                                ButtLR(n)=-1;FixBreak=1;
                                break;
                            end
                        end
                %    end
                end
            end
            if ~RespMade
                [xblah,yblah,buttons] = GetMouse(whichScreen);
                if any(buttons)
                    if buttons(1)==1
                        whichButton=1;
                    elseif buttons(end)==1
                        whichButton=2;
                    end
                    % if par.recordEEG,  io64(eeg.obj,port, par.CD_BUTTONS(whichButton)); WaitSecs(par.eeg_waitTime);  io64(eeg.obj,port,0); end
                    if par.recordEEG,  io64(eeg.obj,port,par.CD_BUTTONS(whichButton)); WaitSecs(par.trigdur);  io64(eeg.obj,port,0); end
                    if par.useEL, Eyelink('Message', ['ButtonPressed']); end
                    ButtT(n) = GetSecs;
                    RespMade=1;
                    nPTBevent = nPTBevent+1;
                    PTBevent(nPTBevent) = par.CD_BUTTONS(whichButton);
                    PTBeventT(nPTBevent) = ButtT(n);
                    % WaitSecs(par.Resp2FBdelay/1000); break;

                end
            end

            %%%%% DISPLAY DOTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            coh = 0;
            if frame >= EvidenceOnsetFrame
                if mod(frame-EvidenceOnsetFrame,par.FPC)<=floor((par.FPC-1)/2) % THIS TINY BIT DOES THE OSCILLATING (ON-OFF) COHERENCE!
                    coh = min(par.cohSet(tgcoh(n)),par.cohSet(tgcoh(n))*(frame-EvidenceOnsetFrame)/numRampFr);
                end
            end

            for df = 1 : dotInfo.numDotField
                % ss is the matrix with 3 sets of dot positions, dots from the last 2 positions and current dot positions
                % Ls picks out the set (e.g., with 5 dots on the screen at a time, 1:5, 6:10, or 11:15)

                % Lthis has the dot positions from 3 frames ago, which is what is then...
                Lthis{df}  = Ls{df}(:,loopi(df));

                % ... moved in the current loop. This is a matrix of random numbers - starting positions of dots not moving coherently.
                this_s{df} = ss{df}(Lthis{df},:);

                % Update the loop pointer
                loopi(df) = loopi(df)+1;

                if loopi(df) == 4,
                    loopi(df) = 1;
                end

                % Make a random selection of dots that will move coherently (0,1 indicator variable for each dot)
                L = rand(ndots(df),1) < coh(df); % this flips a biased coin FOR EACH dot - so variable number of dots move coherently from frame to frame
                % Offset the selected dots
                this_s{df}(L,:) = bsxfun(@plus,this_s{df}(L,:),dxdy{df}(L,:)); % % BSXFUN  Binary Singleton Expansion Function

                if sum(~L) > 0
                    this_s{df}(~L,:) = rand(sum(~L),2);	% get new random locations for the rest
                end

                % Check to see if any positions are greater than 1 or less than 0 which
                % is out of the square aperture, and replace with a dot along one of the
                % edges opposite from the direction of motion.
                N = sum((this_s{df} > 1 | this_s{df} < 0)')' ~= 0; % N indicates for each dot whether out (1) or not (0)
                if sum(N) > 0 % if ANY are out...
                    xdir = sin(pi*dotInfo.dir(df)/180.0);
                    ydir = cos(pi*dotInfo.dir(df)/180.0);
                    % Flip a weighted coin to see which edge to put the replaced dots
                    if rand < abs(xdir)/(abs(xdir) + abs(ydir))
                        this_s{df}(find(N==1),:) = [rand(sum(N),1),(xdir > 0)*ones(sum(N),1)];
                    else
                        this_s{df}(find(N==1),:) = [(ydir < 0)*ones(sum(N),1),rand(sum(N),1)];
                    end
                end

                % Convert to pixel positions for plot 
                this_x{df} = floor(d_ppd(df) * this_s{df});	% pix/ApUnit

                % It assumes that 0 is at the top left, but we want it to be in the
                % center, so shift the dots up and left, which means adding half of the
                % aperture size to both the x and y directions.
                dot_show{df} = (this_x{df} - d_ppd(df)/2)';
            end

            % After all computations, flip to draw dots from the previous loop. For the
            % first time, this doesn't draw anything.

            Screen('FillRect',window, par.rgbFIX, fixRect);
            if frame==1 %%First time dots come on
                if par.diodeflash, Screen('FillRect',window, 255, syncRect); end
                [VBLTimestamp, dotsonT(n)] = Screen('Flip', window,0,dontclear);  %
                if par.recordEEG, io64(eeg.obj,port,par.CD_DOTSON); WaitSecs(par.trigdur);  io64(eeg.obj,port,0); end
                if par.diodeflash, Screen('FillRect',window, par.BGcolor, syncRect); end
                nPTBevent = nPTBevent+1;
                PTBevent(nPTBevent) = par.CD_DOTSON;
                PTBeventT(nPTBevent) = dotsonT(n);
                if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) '_' num2str(par.CD_DOTSON)]); end
            elseif frame==EvidenceOnsetFrame % coherent evidence starts
                if par.diodeflash, Screen('FillRect',window, 255, syncRect); end
                [VBLTimestamp, cohstartT(n)] = Screen('Flip', window,0,dontclear);  %
                if par.recordEEG, io64(eeg.obj,port,par.CD_TG(tgdir(n),tgcoh(n),tgside(n))); WaitSecs(par.trigdur);  io64(eeg.obj,port,0); end
                if par.diodeflash, Screen('FillRect',window, par.BGcolor, syncRect); end
                nPTBevent = nPTBevent+1;
                PTBevent(nPTBevent) = par.CD_TG(tgdir(n),tgcoh(n),tgside(n));
                PTBeventT(nPTBevent) = cohstartT(n);
                if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) 'TG' num2str(par.CD_TG(tgdir(n),tgcoh(n),tgside(n)))]); end
                StimPresented=1;
            else
                Screen('Flip', window,0,dontclear);
                StimPresented=StimPresented + 1;
            end
%                 pause(2); sca; return
            % Now do the actual drawing commands, although nothing is drawn until next
            for df = 1:dotInfo.numDotField
                % NaN out-of-circle dots
                xyDis = dot_show{df};
                outCircle = sqrt(xyDis(1,:).^2 + xyDis(2,:).^2) + dotInfo.dotSize/2 > dotfieldcenter(df,3);
                dots2Display = dot_show{df};
                dots2Display(:,outCircle) = NaN;

                Screen('DrawDots',window,dots2Display,dotInfo.dotSize,dotInfo.dotColor,dotfieldcenter(df,1:2));
                if frame>0, displayed_dots{frame,df} = dots2Display(:,find(~isnan(dots2Display(1,:)))); end
            end

            Screen('DrawingFinished',window,dontclear);   % Tell PTB to get ready while doing computations for next dots presentation

            for df = 1 : dotInfo.numDotField
                % Update the dot position array for the next loop
                ss{df}(Lthis{df}, :) = this_s{df};
            end 

        end
        for df=1:dotInfo.numDotField, for f=1:numFrames, saved_dots{f,df,n} = displayed_dots{f,df}; end, end % save the dots - cause delays ??

    end % end of the trial loop (block)           


    %Save all the clicks during the block:
    RespT = []; RespLR = [];
    [Event, navail] = PsychHID('KbQueueGetEvent' , dev.mouse);
    for k = 1:navail
        [Event, navail2] = PsychHID('KbQueueGetEvent' , dev.mouse);
        if Event.Pressed==1, RespT = [RespT Event.Time]; RespLR = [RespLR Event.Keycode]; end
    end
 
    % Points and Feedback
    I_t = find(PTBevent>=min(par.CD_TG(:)) & PTBevent<=max(par.CD_TG(:))); % find indices of triggers that are target onsets
    targRT = nan(size(I_t)); % leave nan if miss
    targCor = nan(size(I_t)); % correct button?
    totpts = 0; % total points?
    if length(tgdir)~=length(I_t), error('wrong number of targets??'); end
    whichRespAreTargResp = [];
    for n=1:length(I_t)
        nextresp = find(RespT>PTBeventT(I_t(n))+par.minRT,1);
        if ~isempty(nextresp) % IS THERE a response after this target?
            rt = RespT(nextresp) - PTBeventT(I_t(n)); % When was it relative to this target onset?
            if rt < par.maxRT % then it's a legit response to that target
                whichRespAreTargResp = [whichRespAreTargResp nextresp]; % indices of RespT that are target responses
                if RespLR(nextresp)==tgdir(n) % if correct
                    targRT(n) = rt;
                    par.pts(n) = min(par.maxpts,max(0,par.maxpts - (rt-par.ptsDecRT)*par.ptsDecSlope)); 
                    totpts = totpts + par.pts(n);
                    targCor(n) = 1;
                else    % wrong button pressed
                    targCor(n) = 0;
                end
            end
        end
    end
    numMISS = length(find(isnan(targRT)));
    totpts = totpts + par.MISSpenalty*numMISS; % any nans left in targRT must be targets that had no response
    numFA = length(RespT)-length(whichRespAreTargResp);
    totpts = totpts + par.FApenalty*numFA;
    numERR = length(find(targCor == 0));
    totpts = totpts + par.Errpenalty*numERR;
    
 % If it is a practice block and HC, gather the points from LC and make a prediction for LC performance during the task: 
    LCpts =0;
    par.PerfPredLC = 0;
    if par.practice && sum(par.blocktype =='HC')==2 % This is added because Matlab 2015 doesnt have the strings data types yet
            LoCohtr = find(tgcoh==1); % find low Coh trials
            for i=1:length(LoCohtr)
               LCpts = LCpts + par.pts(LoCohtr(i));  % Total Low coherence points in the block
            end
            par.PerfPredLC = LCpts*3; % prediction of average points in a LC 32 trial block with 3 times more LC trials   
    end

              
    % code to check the reward function of RT:
%     r = [0:0.05:2];
%     for i=1:length(r), pts(i) = min(par.maxpts,max(0,par.maxpts - (r(i)-par.ptsDecRT)*par.ptsDecSlope)); end
%     figure; plot(r,pts)

 
    Screen('Flip', window); % make screen blank
    txtstart = 0.1*par.scresw;
    Screen('TextSize', window, 20);
    Screen('DrawText', window, ['On this block you earned a total of ' num2str(round(totpts)) ' points.'], txtstart, 0.10*par.scresh, 255);
    Screen('DrawText', window, ['You clicked the correct button on-time for ' num2str(length(find(~isnan(targRT)))) ' targets '], txtstart, 0.20*par.scresh, 255);
    Screen('DrawText', window, ['out of ' num2str(par.numtrials) '. And you clicked on average ' num2str(mean(targRT(find(~isnan(targRT) & targCor==1)))-par.minRT)], txtstart, 0.25*par.scresh, 255);
    Screen('DrawText', window, ['sec after the target started '], txtstart, 0.3*par.scresh, 255);
    Screen('DrawText', window, ['(the faster, the more points!).'], txtstart, 0.35*par.scresh, 255);
    Screen('DrawText', window, ['And missed ' num2str(numMISS) ' targets (lost 4 points).'], txtstart, 0.40*par.scresh, 255);
    Screen('DrawText', window, ['You clicked the WRONG button for ' num2str(length(find(targRT==0))) ' targets (lost 4 points).'], txtstart, 0.50*par.scresh, 255);
    Screen('DrawText', window, ['And you clicked ' num2str(numFA) ' times IN BETWEEN targets '], txtstart, 0.60*par.scresh, 255);
    Screen('DrawText', window, ['(lost 4 points for each).'], txtstart, 0.65*par.scresh, 255);
    % If doing a HC in a practice block, show predicted performance on LC blocks:
    if  par.practice && sum(par.blocktype =='HC')==2
            Screen('DrawText', window, ['On the HARD trials you earned a total of ' num2str(round(LCpts)) ' points.'], txtstart, 0.75*par.scresh, 255);
            Screen('DrawText', window, ['If you keep this performance, on a hard block' ], txtstart, 0.8*par.scresh, 255);
            Screen('DrawText', window, ['(75% hard trials), you will earn an average of ' num2str(round(par.PerfPredLC)) ' points.' ], txtstart, 0.85*par.scresh, 255);
            Screen('DrawText', window, ['On your next HARD block, try outperform this score!'], txtstart, 0.9*par.scresh, 255);        
    end
            
    Screen('Flip', window);
    
    [clicks,x,y,whichButton] = GetClicks(whichScreen,0);
    
    if par.useEL
        Eyelink('StopRecording');
        Eyelink('CloseFile');
        ELdownloadDataFile
    end
    
    sca;
    ShowCursor
    
    save([par.runID '.mat'],'StartTime','RespT','RespLR','ButtT','saved_dots','tgdir','tgcoh','tgside','PTBevent','PTBeventT','Event','par', 'dotInfo', 'screenInfo','FixonT', 'dotsonT', 'cohstartT')
    save parCDots par

catch
    sca;
    ShowCursor
    ple
    
    if par.useEL
        Eyelink('StopRecording');
        Eyelink('CloseFile');
    end
   
    %     ListenChar(0);
end

try
    PsychHID('KbQueueFlush', dev.mouse);
    PsychHID('KbQueueStop', dev.mouse);
    PsychHID('KbQueueRelease', dev.mouse);
catch
end
