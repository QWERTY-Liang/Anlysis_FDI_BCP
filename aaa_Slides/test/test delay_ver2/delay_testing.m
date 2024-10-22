%%RealTime Data Streaming with Delsys SDK

% Copyright (C) 2011 Delsys, Inc.
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, and distribute the
% Software, and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.

function Training_MVC
clear all
close all

clc

try
    
    Screen('Preference', 'SkipSyncTests', 1);   %%This should be zero for main experiment %%%%%%%%%%%%%%%%change back to zero
    Screen('CloseAll')
    
    %%Initial parameters for dialog box
    par.pID='P'; par.MVC=0; par.numtrials=20; par.thresh= [0.1];   % on first run of task, set defaults to put in dialog box
    
    dlg_title = 'Exp Parameters';
    while 1
        
        prompt = {'Enter PARTICIPANT', 'RUN/TASK IDENTIFIER:','MVC? (1=max, 0=no, 2= rest)', ['Number of trials?'], ['Threshold? % of MVC']}; %
        def = {par.pID, '', num2str(par.MVC), num2str(par.numtrials),num2str(par.thresh)};
        answer = inputdlg(prompt,dlg_title,1,def);
        par.pID = answer{1};
        par.run = answer{2};
        par.runID = [par.pID '_' par.run];
        par.MVC = str2num(answer{3});
        par.numtrials = str2num(answer{4});
        par.thresh = str2num(answer{5});
        
        
        if par.MVC==1 && exist([par.pID '_MVC.mat'],'file'),
            dlg_title = [par.pID '_MVC.mat EXISTS ALREADY - CHOOSE ANOTHER, OR DELETE THAT ONE IF IT IS RUBBISH']
        elseif par.MVC==2 && exist([par.pID '_Rest.mat'],'file'),
            dlg_title = [par.pID '_Rest.mat EXISTS ALREADY - CHOOSE ANOTHER, OR DELETE THAT ONE IF IT IS RUBBISH']
        elseif ~ par.MVC && exist([par.runID '_train.mat'],'file'),
            dlg_title = [par.runID '_train.mat EXISTS ALREADY - CHOOSE ANOTHER, OR DELETE THAT ONE IF IT IS RUBBISH']
        else
            break
        end
        
        
    end
    
    percentage_of_MVC = par.thresh;
    par.ITI=1;
    if ~par.MVC
        [MVC1,MVC2] = readMVC([par.pID '_MVC.mat']);
        [Rest1,Rest2] = readRest([par.pID '_Rest.mat']);
    end
    %%Important task parameters
    windowtime =250e-3;  %% in seconds the length of EMG window for MAV calculation
    EMGWindowLength = ceil(1926*windowtime); %%In sample points.
    maxTime = 3;  %%maximum time for trial in seconds
    time_over_th = 4; %Number of refreshes to requre over threshold
    par.videoFrate = 85;
    %%EMG filter
    w1 = 10/(1926/2);
    w2 = 500/(1926/2);
    [bcoef, acoef] = butter(4, [w1 w2]);
    
    
    if par.MVC
        maxTime = 1;
    end
    
    
    %% Initialize Psychtoolbox
    
    PsychDefaultSetup(2);
    screens = Screen('Screens');
    screenNumber = max(screens);
   
   
    % Open window
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, [0 0 0]);
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
    
    par.hz = Screen('FrameRate', 0, 1); % find out the refresh rate as well
    disp(['Screen Refresh Rate ' num2str(par.hz)]);
        % Check if screen parameters are correct
%     if scres~= [1920,1080] %[1280 720] %[1024,768] %[1280 720]this aspect ratio is according to TUF monitor % %[1920 1080]   %[2736,1824] % [1280,1024]%  DESIRED SR = 1024, 768. Changed to 1920 to run on LS_NOVA
%         error(['The monitor is NOT SET to the Screen Resolution. Please change it.'])
%     end
    
    if abs(par.videoFrate-par.hz)>1
        error(['The monitor is NOT SET to the desired frame rate of ' num2str(par.videoFrate) ' Hz. Change it.'])
    end
    % Define rectangle dimensions
    width = 200;
    spacing = 475;
    
    % Calculate rectangle positions
    centerX = round(0.5 * RectWidth(Screen('Rect', screenNumber)));
    rectLeftX = centerX - spacing - round(0.5 * width);
    rectRightX = centerX + spacing - round(0.5 * width);
    rectBottomY = RectHeight(Screen('Rect', screenNumber));
    
    
    % Calculate threshold line position
    
    thresholdY = (1-percentage_of_MVC) * RectHeight(Screen('Rect', screenNumber)); %(bottem of screen = 100%!!!!!!!!!)
    thresholdX = RectWidth(Screen('Rect', screenNumber));
    
    %% Load starting screen
    
    handoptions = {'left', 'right'};
    hand = randperm(2); %will randomly order the values in handoptions (l, r) or (r,l)
    hand = hand(1);  %assign hand the first cell (l or r)
    
    %creast intro to be shown on screen
%     intro1 = 'Trial Commencing';
%     if par.MVC
%         intro1 = 'MVC Trial Commencing';
%     end
%     intro2 = '\n Relax Muscles';
%     intro3 = sprintf(['\n Trial Number: 1 - ' handoptions{hand}]);
%     
%     Screen('TextSize', window, 100);
%     DrawFormattedText(window,[intro1 intro2 intro3] ,'center', screenYpixels * 0.5, [1 1 1]); %centre words on screen in white
%     
    % Flip to the screen
    Screen('Flip', window);
    
    tockp=[];
    tic
    for i=1:100
        
        Screen('flip',window)
        tockp= [ tockp; toc];
    end
    if abs(1/(median(diff(tockp))) - par.hz)>1, error('refresh rate off - try restarting matlab or the computer'), end
  
    
    %wait 5s
    WaitSecs(1); %change to 5
    
    %%
    %%EMG setup
    % CHANGE THIS TO THE IP OF THE COMPUTER RUNNING THE TRIGNO CONTROL UTILITY
    HOST_IP = '137.43.149.175';  %'192.168.1.60'
    par.monitorwidth_cm = 60;
    load('gammafn_UCDRight_3par.mat');


    %%%%test code delay#
    trig.ioObj = io64;
    %initialize the inpoutx64 system driver
    status = io64(trig.ioObj);
    if(status == 0)
        disp("EMG triggers ready")
    end
    trig.address = hex2dec('4FF8');
    
    io64(trig.ioObj, trig.address, 0);
        %% EMG setup
    % CHANGE THIS TO THE IP OF THE COMPUTER RUNNING THE TRIGNO CONTROL UTILITY
    % if booth_name == 'L'   %%%Simon's original setup
    %     HOST_IP = '137.43.148.180';    %%Needed for real-time EMG recording
    %     par.monitorwidth_cm = 54.3;
    %     load('gammafn_UCDnew_3par.mat');
    % elseif booth_name == 'R'  %%New booth setup
    %     HOST_IP = '137.43.149.175';    %%Needed for real-time EMG recording
    %     par.monitorwidth_cm = 60;
    %     load('gammafn_UCDRight_3par.mat');
    % end

    
    %% Create the required objects
    
    
    
    %TCPIP Connection to stream EMG Data
    interfaceObjectEMG = tcpip(HOST_IP,50041); %original reg trignos port=41, need port=43 for dual trignos
    interfaceObjectEMG.InputBufferSize = 6400*18;   %600 samples
    num_buffer_samples = interfaceObjectEMG.InputBufferSize/64;
    %TCPIP Connection to communicate with SDK, send/receive commands
    commObject = tcpip(HOST_IP,50040);
    
    
    
    %% Open the COM interface, determine RATE
    
    fopen(commObject);
    
    pause(1);
    fread(commObject,commObject.BytesAvailable);
    fprintf(commObject, sprintf(['RATE 2000\r\n\r']));
    pause(1);
    fread(commObject,commObject.BytesAvailable);
    fprintf(commObject, sprintf(['RATE?\r\n\r']));
    pause(1);
    data = fread(commObject,commObject.BytesAvailable);
    
    emgRate = strtrim(char(data'));
    % if(strcmp(emgRate, '1925.926'))
    %     rateAdjustedEmgBytesToRead=1664;
    % else
    %     rateAdjustedEmgBytesToRead=1728;
    % end
    
    rateAdjustedEmgBytesToRead = 640;  %%will read in either 20 or 30 samples with the 85Hz refresh rate
    %This just needs to be a multiple of 64, since we are not using the
    %accelerometer. It's the minimum number of bytes to read. 64 bytes is
    %one sample.
    
    
    %drawnow
    
    
    %pause(1);
    %%
    % Open the interface object
    try
        fopen(interfaceObjectEMG);
    catch
        sca
        Screen('CloseAll');
        keyboard
        
        error('CONNECTION ERROR: Please start the Delsys Trigno Control Application and try again');
    end
    
    
    
    %%
    % Send the commands to start data streaming
    WaitSecs(par.ITI);
    
    fprintf(commObject, sprintf(['START\r\n\r']));
    
    WaitSecs(par.ITI);
    pause(2);
    
    outcome  = nan(1,par.numtrials); % outcome = 0 for timeout, 1 for successful left hand trial, 2 for successful right hand trial
    start_time_emg  = nan(1,par.numtrials);
   pre_buffer_time = nan(1,par.numtrials);
    if par.MVC
        MVC1 = nan(1,par.numtrials);
        MVC2 = nan(1,par.numtrials);
        Rest1 = nan(1,par.numtrials);
        Rest2 = nan(1,par.numtrials);
    end
    cue_hand = nan(1,par.numtrials);
    trialstartT = nan(1,par.numtrials);
    trialRT = nan(1,par.numtrials);
    EMG_buffer=[];
    buffer_size =1926; %%make this a par?
    %%%%%%%%%%
      Delay = [];
    %%%%%%%%%%
    for trialnum = 1: par.numtrials

        io64(trig.ioObj, trig.address, 0);
        bytesAvailRec  = [];
        data_arrayBaseline=[];
        data_arrayEMG=[];
        baselineread{trialnum}=[];
        baselinesamples{trialnum}=[];
        num_samples{trialnum}=nan(1, ceil(par.hz*maxTime));
        MAV{trialnum}=nan(ceil(par.hz*maxTime),2);
        beforeBuffRead_time{trialnum}=nan(1, ceil(par.hz*maxTime));
        afterBuffRead_time{trialnum}=nan(1, ceil(par.hz*maxTime));
        bytesAvailRec{trialnum}=nan(1, ceil(par.hz*maxTime));
        %counts for determining how much time elapses while rect is over the threshold
        count1=0;
        count2=0;
        EMGlevel1 = 0;
        EMGlevel2 = 0;
        cue_hand(trialnum) = hand;
        
        %%Let's read in baseline EMG for 1 sec or so during instructions to make sure buffer
        %%cleared
        
        pre_buffer_time(trialnum) = GetSecs;
        keepgoing=1;
        firstframe = 1;
        if trialnum==1, baselinetime = 6; else baselinetime = 3; end;
        while keepgoing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if firstframe
                t1 = GetSecs();
                io64(trig.ioObj, trig.address, 1);
                while 1
                    bytesReadyTMP = interfaceObjectEMG.BytesAvailable;
                    %bytesAvailRec = [bytesAvailRec bytesReady];
                    bytesReadyTMP = bytesReadyTMP - mod(bytesReadyTMP, rateAdjustedEmgBytesToRead);%%1664
                    if ~(bytesReadyTMP == 0)
                        % beforeBuffRead_time{n}((y-1)*numRefreshPerCycle + i) = GetSecs;
                        TMP  = cast(fread(interfaceObjectEMG,bytesReadyTMP), 'uint8');
                        baselineread{trialnum} = [ baselineread{trialnum}; GetSecs];
                        TMP = typecast(TMP, 'single');
                        data_ch15 = TMP(2:16:end);
                        
                        baselinesamples{trialnum} = [ baselinesamples{trialnum}; length(data_ch15)];

                    else
                         baselinesamples{trialnum} = [ baselinesamples{trialnum}; 0];
                    end
                    if min(data_ch15) < -1
                        t2 = GetSecs();
                        Delay(end+1) = t2 - t1;
                        io64(trig.ioObj, trig.address, 0);
                        break
                    end
    
                end
                firstframe = 0;
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            bytesReadyTMP = interfaceObjectEMG.BytesAvailable;
            %bytesAvailRec = [bytesAvailRec bytesReady];
            bytesReadyTMP = bytesReadyTMP - mod(bytesReadyTMP, rateAdjustedEmgBytesToRead);%%1664
            
            if ~(bytesReadyTMP == 0)
                % beforeBuffRead_time{n}((y-1)*numRefreshPerCycle + i) = GetSecs;
                TMP  = cast(fread(interfaceObjectEMG,bytesReadyTMP), 'uint8');
                baselineread{trialnum} = [ baselineread{trialnum}; GetSecs];
                TMP = typecast(TMP, 'single');
                data_ch1 = TMP(11:16:end);
                data_ch3 = TMP(2:16:end);
                data_arrayBaseline = [data_arrayBaseline; [data_ch1,data_ch3]];
                EMG_buffer = [EMG_buffer; [data_ch1,data_ch3]];
                
                baselinesamples{trialnum} = [ baselinesamples{trialnum}; length(data_ch1)];
            else
                 baselinesamples{trialnum} = [ baselinesamples{trialnum}; 0];
            end
            %show load screen
            line1 = 'New Trial Commencing';
            if par.MVC == 1
                line1 = 'New MVC Trial Commencing';
            elseif par.MVC == 2
                line1 = 'New Rest Trial Commencing';
            end
            line2 = '\n Relax Muscles';
            line3 = sprintf(['\n Trial Number: %d - ' handoptions{hand}],trialnum);
            Screen('TextSize', window, 100);
            DrawFormattedText(window, [line1 line2 line3],'center', screenYpixels * 0.5, [1 1 1]);
            
            % 
             if GetSecs-pre_buffer_time(trialnum) > 10
                error('EMG problem?')

            end
     


            if size(EMG_buffer,1) > buffer_size
                if EMG_buffer(end-buffer_size,1) ~= 0 && GetSecs-pre_buffer_time(trialnum) > 3 && length(data_ch1) < num_buffer_samples
                    keepgoing = 0;
                end
            end
            
            % Flip to the screen
            Screen('Flip', window);
        end
        start_time_emg(trialnum) = GetSecs;

        %io64(trig.ioObj, trig.address, 1);%%%%%%%%%%%%%%%%%%%%%%%%

        EMG_buffer = double(EMG_buffer(end-(buffer_size-1):end,:));
        initial_EMG_buffer = EMG_buffer;


        for ind = 1:ceil(par.hz*maxTime)

            % if ind>par.hz*1.5
            % 
            %     io64(trig.ioObj, trig.address, 1);%%%%%%%%%%%%%%%
            % else
            % 
            %     io64(trig.ioObj, trig.address, 0);%
            % end



            if  ((count1>= time_over_th) & hand==1) | ((count2 >= time_over_th)  & hand==2) %if rect is over line > 400 msec
                
                Screen('TextSize', window, 30);
                DrawFormattedText(window, 'Well done','center', screenYpixels * 0.5, [1 1 1]);
                Screen('Flip', window);
                if isnan(outcome(trialnum))
                    
                    outcome(trialnum) = hand;
                end
                
                break
            else
                
                bytesReady = interfaceObjectEMG.BytesAvailable;
                bytesAvailRec{trialnum}(ind) = bytesReady;
                bytesReady = bytesReady - mod(bytesReady, rateAdjustedEmgBytesToRead);%%1664
                if ~(bytesReady == 0)
                    
                    beforeBuffRead_time{trialnum}(ind) = GetSecs;
                    data = cast(fread(interfaceObjectEMG,bytesReady), 'uint8');
                    data = typecast(data, 'single');
                    afterBuffRead_time{trialnum}( ind) = GetSecs;
                    %bytesReady = interfaceObjectEMG.BytesAvailable;
                    %bytesAvailRec2 = [bytesAvailRec2 bytesReady];
                    
                    data_ch1 = data(11:16:end); %left hand - will be plotted as rect1
                    data_ch3 = data(2:16:end); %right hand (using channel3 - will be plotted as rect2
                     num_samples{trialnum}(ind) = length(data_ch1);
                    data_arrayEMG = [data_arrayEMG; [data_ch1,data_ch3]];
                    EMG_buffer = [EMG_buffer(length(data_ch1)+1:end,:);  double([data_ch1, data_ch3])];
                end
             %   if length(data_arrayEMG)>=EMGWindowLength 
             
                    
%                     %filter
                    filteredsignal1 = filtfilt(bcoef, acoef, (EMG_buffer(:,1)));
                    filteredsignal2 = filtfilt(bcoef, acoef, (EMG_buffer(:,2)));
                   
                    %    EMGwindow{n}(emgind, :, :) = [filteredsignal1 filteredsignal2]; %Probably to big to save
                    %these are the values we want to plot as rect heights
                    EMGlevel1 = mean(abs(filteredsignal1((end-EMGWindowLength+1):end))); %left hand/rect
                    EMGlevel2 = mean(abs(filteredsignal2((end-EMGWindowLength+1):end))); %right hand/rect



                    MAV{trialnum}(ind,:) = [EMGlevel1 EMGlevel2];
                    
                    
                    
                    if ~par.MVC
                        % Generate new random heights for rectangles
                        % Normalise EMG
                        rectHeight1 = max([((EMGlevel1-Rest1)/(MVC1-Rest1))*screenYpixels, 0]);
                        rectHeight2 = max([((EMGlevel2-Rest2)/(MVC2-Rest2))*screenYpixels, 0]);
                        
                        %%if over threshold
                        if rectBottomY-rectHeight1<=thresholdY
                            colour1 =  double([0.2 0.9 0.2]); %green
                            count1=count1+1;
                            
                            if hand==1  && isnan(trialRT(trialnum))
                                trialRT(trialnum)= GetSecs;
                            end
                        else
                            colour1 = double([1 0 0.1]); %red
                            count1=0;
                            if hand==1
                                trialRT(trialnum) = nan;
                            end
                        end
                        
                        if rectBottomY-rectHeight2<=thresholdY
                            colour2 =  double([0.2 0.9 0.2]); %green
                            count2=count2+1;
                            if hand==2 && isnan(trialRT(trialnum))
                                trialRT(trialnum)= GetSecs;
                            end
                        else
                            colour2 = double([1 0 0.1]);  %red
                            count2=0;
                            if hand==2
                                trialRT(trialnum) = nan;
                            end
                        end
                        %
                        
                        % Draw rectangles
                        rect1 = double([rectLeftX, rectBottomY - rectHeight1, rectLeftX + width, rectBottomY]);
                        rect2 = double([rectRightX, rectBottomY - rectHeight2, rectRightX + width, rectBottomY]);
                        
                        Screen('FillRect', window, colour1, rect1);
                        Screen('FillRect', window, colour2, rect2);
                        
                        
                        % Label rectangles and counters
                        label1 = sprintf('left hand ');
                        label2 = sprintf('right hand');
                        Screen('TextSize', window, 30);
                        DrawFormattedText(window, label1, rectLeftX-200, rectBottomY -30, [1 1 1]);
                        DrawFormattedText(window, label2, rectRightX+250, rectBottomY - 30, [1 1 1]);
                        
                        % Draw threshold line and label
                        Screen('DrawLine', window, [1 1 1], 0, thresholdY, thresholdX, thresholdY,10);
                        DrawFormattedText(window, 'threshold', thresholdX - 200, thresholdY - 30, [1 1 1]);
                        
                        %Screen('Flip', window,0.9);
                        if isnan(trialstartT(trialnum))
                            trialstartT(trialnum) =  Screen('Flip', window);  %%Trial startT is when the interface is visible to the participant. start_time_emg is when the trial EMG starts getting saved
                       
                        else
                         Screen('Flip', window);
                        end
                        
                    else
                        Screen('Flip', window);
                    end
                %else
                     Screen('Flip', window);
                
            end
        end
        
        
        
        
        
        
        WaitSecs(0.5)    %freeze bar for 2ss for recognition of success
        
        
        RT(trialnum) = trialRT(trialnum) - trialstartT(trialnum);
        
        if trialnum<par.numtrials  %%instructions for next trial
            hand = 3-hand; %switch hand
            
            
            
            %wait 5s
            
        end
        
        %do we want to use 80/90% max 'quantile' func
        if par.MVC == 1
            MVC1(trialnum) = max( MAV{trialnum}(:,1));
            MVC2(trialnum) = max( MAV{trialnum}(:,2));
            
        elseif par.MVC == 2
            Rest1(trialnum) = nanmedian( MAV{trialnum}(:,1));
            Rest2(trialnum) = nanmedian( MAV{trialnum}(:,2));
            
        end
        
        FullEMG{trialnum} = [initial_EMG_buffer; data_arrayEMG];
        BaselineEMG{trialnum} = data_arrayBaseline;
    end
    
    %save important variables before we fully close the interface
    if par.MVC ==1
        MVC1 = max(MVC1);
        MVC2 = max(MVC2);
        save([par.pID '_MVC'], 'Delay', 'cue_hand', 'MAV', 'FullEMG', 'BaselineEMG', 'trialstartT', 'EMGWindowLength', 'maxTime', 'MVC1', 'MVC2', 'baselineread', 'start_time_emg','afterBuffRead_time','beforeBuffRead_time','start_time_emg', 'pre_buffer_time');
    elseif par.MVC ==2
        Rest1 = min(Rest1);
        Rest2 = min(Rest2);
        save([par.pID '_Rest'], 'cue_hand', 'MAV', 'FullEMG', 'BaselineEMG', 'trialstartT', 'EMGWindowLength', 'maxTime', 'Rest1', 'Rest2', 'baselineread', 'start_time_emg','afterBuffRead_time','beforeBuffRead_time','start_time_emg', 'pre_buffer_time');
    else
        save([par.runID '_train.mat'], 'par', 'cue_hand', 'outcome', 'MAV', 'FullEMG', 'BaselineEMG', 'trialstartT', 'trialRT', 'RT', 'EMGWindowLength', 'maxTime', 'time_over_th', 'MVC1', 'MVC2', 'percentage_of_MVC', 'num_samples', 'beforeBuffRead_time', 'afterBuffRead_time', 'baselineread', 'baselinesamples', 'start_time_emg', 'pre_buffer_time' , 'Rest1', 'Rest2', 'bytesAvailRec');
    end
    Screen('CloseAll');
    
    
    
catch
    sca
    Screen('CloseAll');
    ple
    keyboard
    
    ShowCursor
    
    
    % Clean up the network objects
    if isvalid(interfaceObjectEMG)
        fclose(interfaceObjectEMG);
        delete(interfaceObjectEMG);
        clear interfaceObjectEMG;
    end
    
    
    
    if isvalid(commObject);
        fclose(commObject);
        delete(commObject);
        clear commObject;
    end
end


% Clean up the network objects
if isvalid(interfaceObjectEMG)
    fclose(interfaceObjectEMG);
    delete(interfaceObjectEMG);
    clear interfaceObjectEMG;
end



if isvalid(commObject);
    fclose(commObject);
    delete(commObject);
    clear commObject;
end
keyboard
