%%
%这个EOL与RL Sl都能做
%1.每人单独分析
%2. 加了'调整allbehaviour'以及检查
%% Toolbox requirements:
clc
clear
addpath('G:\My Drive\Phd\EEGLAB\eeglab-develop');% EEGlab toolbox
addpath(genpath('G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2'));% current folder

%% Set experimental analysis parameters
exp.sub_id = [1,2,3,4,5,6];
%exp.sub_id = [1];
[exp] = TLBEM1_setup(exp);

%% Run analysis
eeglab

%还需手动变换分类条件！！下次更新

for fband=1:5 %控制不同的频率段组合
    %ff = find((F>8 & F<21.1) | (F>21.2 & F<30)); % the indices of F that cover the spectral range/band of interest. Let's say we're interested in Mu and Beta bands combined. Note here I'm avoiding the SSVEP frequency (18.75hz in this example)
    %ff = find((F>13 & F<21.1) | (F>21.2 & F<30)); %只有Beta
    %ff = find(F>8 & F<13 ); %只有Mu
    %ff = find(F>1 & F<4 ); % Delta analysis
    %ff = find(F>4 & F<7 );
    disp(['now is: ' num2str(fband) 'frequency band'] );
    for e=1:3 % 2 是RL;3是EOL\ 1 是SL
        disp(['now is: ' num2str(e) 'epoch type'] );
        for sub=1:6% 每个人分别看
            epoch = exp.epochs{e}; %选择切分方法

            %% 1. Load rl_BB pre-cue-188
            % col 1: subject num
            % col 2: block num
            % col 3: contrast (high or low contrast)
            % col 4: muscle (4=FDI, 44=BCP)
            % col 5: trial outcome (1=correct, 2=error, 4=no response, 3 too early, 6 slow, 5 wrong muscle)
            % col 6: RT in sec
            % col 7: evshowtime =[]   % 0.2sec additional time after response
            % col 8: participant response, 1 = left, 2 = right 0= on response
            % col 9: correct response, 1 = left, 2 = right
            % col 10: first tilt used for SSVEP
            %证据前0.188秒
            % cue前0.188秒 col:12
            % 整段平均 col:13
            % response 前0.188秒 col:14
            addpath(genpath(exp.finalpath));

            %加载行为数据
            load([exp.behpath exp.name 'EMG_brusts_update'])
            AllBehaviour=AllBehaviour_EMGonsite((sub-1)*1024+1:(sub)*1024,:);
            clear AllBehaviour_EMGonsite


            %测试无baseline情况
            %EEG = pop_loadset([exp.finalpath 'csd_acICAri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);
            %测试有baseline情况 (5 1024)
            EEG = pop_loadset([exp.finalpath 'csd_abb_cICAri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.set']);

            filename=['Rejected_bb_cICAri' epoch '_' exp.filterLab 'aac' exp.name num2str(sub) '.mat'];
            load(filename)

            if length(bTrial_ind)-length(bTrial_num)~=EEG.trials
                disp(sub);
                error('document not match')
            else
                disp('length checked')
            end

            %% 调整 ALLbehaviour
            for ii=flip(1:length(bTrial_ind))
                if bTrial_ind(ii)==1
                    AllBehaviour(ii,:)=[];
                end
            end

            %% STFT prepare
            fs=EEG.srate;%512
            t=EEG.times;
            erp=EEG.data;

            % epoch_limits_msTG = [-188*4 188*12];    % Target-locked epoch. We could use the same windows as were used for artifact rejection, or we might want different windows here - if so, don't forget that artifacts were only checked for in the window specified in IdentifyBadTrials!
            % tts = round(epoch_limits_msTG(1)/1000*fs):round(epoch_limits_msTG(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
            % tt = tts*1000/fs; % hence timebase in milliseconds, for plotting etc
            % epoch_limits_msR =  [-188*3 188];    % Response-locked epoch; again using integer number of SSVEP cycles (1000/18.75 = 53.33)
            % trs = round(epoch_limits_msR(1)/1000*fs):round(epoch_limits_msR(2)/1000*fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
            % tr = trs*1000/fs; % hence timebase in milliseconds


            %这里可以分别看mu和 beta
            % Let's say we also want to compute spectral anmplitude as a function of time in the Mu/beta bands (reflects motor preparation when measured
            % over motor cortices). One way to do such time-frequency analysis is the short-Time Fourier Transform - taking FFTs in a sliding window
            % across time. First define parameters of this analysis:
            if e==2 || e==3
                fftlen = round(fs/21.25*6); % Window of how many sample points? If there is an SSVEP involved, whether or not you are interested in analyzing it, it is good to have all power related to the SSVEP isolated in a single frequency bin. This happens when you choose a window length that is an integer number of SSVEP cycles.
            elseif e==1
                fftlen = round(fs/21.25*12);
            end


            F = [0:fftlen-1]*fs/fftlen; % frequency scale, given window length (remember resolution = 1/window-duration)
            if fband==1
                ff = find((F>8 & F<21.1) | (F>21.2 & F<30)); % the indices of F that cover the spectral range/band of interest. Let's say we're interested in Mu and Beta bands combined. Note here I'm avoiding the SSVEP frequency (18.75hz in this example)
                fbandname='MuBeta';
            elseif fband==2
                ff = find((F>13 & F<21.1) | (F>21.2 & F<30)); %只有Beta
                fbandname='Beta';
            elseif fband==3
                ff = find(F>8 & F<13 ); %只有Mu
                fbandname='Mu';
            elseif fband==4
                ff = find(F>1 & F<4 ); % Delta analysis
                fbandname='Delta';
            elseif fband==5
                ff = find(F>4 & F<7 );
                fbandname='Theta';

            end

            if e==2 || e==3
                Ts = [-2700:20:800];%[-188*4:47:188*12]; % in msec, centered on what times do you want to measure spectral amplitude? i.e., where to center each consecutive window in time
                % Tr = [-188*3:47:188]; % for response-locked
            elseif e==1
                Ts = [-950:20:2700];
            end
            % Now in the rest of the code in the loop, we turn to computing and averaging Mu/Beta amplitude
            % compute short time fourier transform (STFT) for each single trial:
            STFT = []; % initialise
            % for the stimulus-locked STFT, we can compute ffts across all trials at once - that's why there is no 'for n=1:ntr' loop.
            for tt=1:length(Ts) % for each STFT timepoint (window centre) we'll compute the FFT for all trials at once
                [blah,samp] = min(abs(t-Ts(tt))); % find the sample point in the ERP epoch corresponding to the centre of the current FFT window
                spec = abs(fft(erp(:,samp-round(fftlen/2)+[1:fftlen],:),[],2))./(fftlen/2); % compute the magnitude of FFT (and scale it so that it's amplitude in same units as EEG
                % Save the result for each trial, just like the matrix 'erp'
                STFT(:,tt,:) = mean(spec(:,ff,:),2);
                disp(tt)
            end

            %% 5.1 sorting left vs right (Correct vs wrong)
            if e==2
                RT= AllBehaviour(:,6); % 如果RL则RT分类

            elseif e==3
                RT= AllBehaviour(:,11); % 如果EOL则EMG onsite 分类
                for j=1:length(RT)
                    if RT(j)>3
                        RT(j)=NaN;
                    end

                end
            elseif e==1

                RT= AllBehaviour(:,6); % 如果SL则 end of trial 分类%%%% 这里有bug, 改为RT分类

            end


            RT33=prctile(RT,100/3);%三等分低点
            RT66=prctile(RT,200/3);%三等分高点
            % left_channel = 115;  % C3通道号
            % right_channel = 54;
            selected_trials_1 = 999*ones(length(AllBehaviour),3);  % left C4-C3
            selected_trials_2 = 999*ones(length(AllBehaviour),3);  % right C3-C4
            selected_trials_11 = zeros(length(AllBehaviour),3);  % left C4-C3
            selected_trials_22 = zeros(length(AllBehaviour),3);  % right C3-C4

            for RTsplit=1:3
                % selected_trials_1 = 999*ones(length(AllBehaviour),RTsplit);  % left C4-C3
                % selected_trials_2 = 999*ones(length(AllBehaviour),RTsplit);  % right C3-C4
                count=0;
                % switch%%%%%%%%%%calculate correct first
                %     case 'CW_SL_b' %correct vs wrong; pre cue SL
                for i=1:length(AllBehaviour)
                    %select trials
                    if AllBehaviour(i,8)==1&& AllBehaviour(i,9)==2% &&AllBehaviour(i,4)==4% &&&& AllBehaviour(i,3)==0.14%&& AllBehaviour(i,9)==1%&& AllBehaviour(i,3)==0.14%&& AllBehaviour(i,4)==4% && AllBehaviour(i,4)==44
                        selected_trials_1(i,RTsplit)=1;
                    elseif AllBehaviour(i,8)==2&& AllBehaviour(i,9)==1%&&AllBehaviour(i,4)==4%&& AllBehaviour(i,3)==0.14 %&& AllBehaviour(i,9)==2%&& AllBehaviour(i,3)==0.14%&& AllBehaviour(i,4)==4%&& AllBehaviour(i,4)==44
                        selected_trials_2(i,RTsplit)=1;
                    end
                    %exclude invalide(too early or wrong muscle)
                    if AllBehaviour(i,5)==3 || AllBehaviour(i,5)==5|| AllBehaviour(i,5)==4%|| AllBehaviour(i,5)==6
                        selected_trials_1(i,RTsplit)=0;%change back to unselected
                        selected_trials_2(i,RTsplit)=0;
                        count=count+1; %计数丢掉的trial
                    end

                    %把不属于时间的剔除
                    if RTsplit==1
                        if RT(i)>=RT33%只保留偏小三等分点下的
                            selected_trials_1(i,RTsplit)=0;%change back to unselected
                            selected_trials_2(i,RTsplit)=0;
                        end
                        trl1_1=find(selected_trials_1(:,RTsplit)==1);
                        trl2_1=find(selected_trials_2(:,RTsplit)==1);

                    elseif RTsplit==2
                        if RT(i)<=RT33 || RT(i)>=RT66%只保留两个三等分点中的
                            selected_trials_1(i,RTsplit)=0;%change back to unselected
                            selected_trials_2(i,RTsplit)=0;
                        end

                        trl1_2=find(selected_trials_1(:,RTsplit)==1);
                        trl2_2=find(selected_trials_2(:,RTsplit)==1);
                    elseif RTsplit==3
                        if RT(i)<=RT66%只保留偏da三等分点上的
                            selected_trials_1(i,RTsplit)=0;%change back to unselected
                            selected_trials_2(i,RTsplit)=0;
                        end
                        trl1_3=find(selected_trials_1(:,RTsplit)==1);
                        trl2_3=find(selected_trials_2(:,RTsplit)==1);

                    end

                end



                %low
                % selected_trials_11 = zeros(length(AllBehaviour),3);  % left C4-C3
                % selected_trials_22 = zeros(length(AllBehaviour),3);  % right C3-C4
                count=0;


                % switch%%%%%%%%%%calculate correct first
                %     case 'CW_SL_b' %correct vs wrong; pre cue SL
                for i=1:length(AllBehaviour)
                    %select trials
                    if AllBehaviour(i,8)==1&& AllBehaviour(i,9)==2% &&AllBehaviour(i,4)==44%4&&AllBehaviour(i,3)==0.07 %&& AllBehaviour(i,9)==2%&&AllBehaviour(i,3)==0.07% AllBehaviour(i,4)==44%&& AllBehaviour(i,3)==0.07
                        selected_trials_11(i,RTsplit)=1;
                    elseif AllBehaviour(i,8)==2&& AllBehaviour(i,9)==1% &&AllBehaviour(i,4)==44%&&AllBehaviour(i,3)==0.07 %&& AllBehaviour(i,9)==1%&& AllBehaviour(i,3)==0.07%AllBehaviour(i,4)==44%&& AllBehaviour(i,3)==0.07
                        selected_trials_22(i,RTsplit)=1;
                    end
                    %exclude invalide(too early or wrong muscle)
                    if AllBehaviour(i,5)==3 || AllBehaviour(i,5)==5|| AllBehaviour(i,5)==4
                        selected_trials_11(i,RTsplit)=0;%change back to unselected
                        selected_trials_22(i,RTsplit)=0;
                        count=count+1;
                    end

                    %把不属于时间的剔除
                    if RTsplit==1
                        if RT(i)>=RT33%只保留偏小三等分点下的
                            selected_trials_11(i,RTsplit)=0;%change back to unselected
                            selected_trials_22(i,RTsplit)=0;
                        end

                        trl11_1=find(selected_trials_11(:,RTsplit)==1);
                        trl22_1=find(selected_trials_22(:,RTsplit)==1);


                    elseif RTsplit==2
                        if RT(i)<=RT33 || RT(i)>=RT66%只保留两个三等分点中的
                            selected_trials_11(i,RTsplit)=0;%change back to unselected
                            selected_trials_22(i,RTsplit)=0;
                        end


                        trl11_2=find(selected_trials_11(:,RTsplit)==1);
                        trl22_2=find(selected_trials_22(:,RTsplit)==1);
                    elseif RTsplit==3
                        if RT(i)<=RT66%只保留偏da三等分点上的
                            selected_trials_11(i,RTsplit)=0;%change back to unselected
                            selected_trials_22(i,RTsplit)=0;
                        end

                        trl11_3=find(selected_trials_11(:,RTsplit)==1);
                        trl22_3=find(selected_trials_22(:,RTsplit)==1);

                    end

                end





            end
            %% compute
            %
            avMB1_1(:,:,sub) = mean(STFT(:,:,trl1_1(:)),3);
            avMB2_1(:,:,sub) = mean(STFT(:,:,trl2_1(:)),3);
            avMB11_1(:,:,sub) = mean(STFT(:,:,trl11_1(:)),3);% wrong-left
            avMB22_1(:,:,sub) = mean(STFT(:,:,trl22_1(:)),3);% wrong-right

            avMB1_2(:,:,sub) = mean(STFT(:,:,trl1_2(:)),3);
            avMB2_2(:,:,sub) = mean(STFT(:,:,trl2_2(:)),3);
            avMB11_2(:,:,sub) = mean(STFT(:,:,trl11_2(:)),3);% wrong-left
            avMB22_2(:,:,sub) = mean(STFT(:,:,trl22_2(:)),3);% wrong-right

            avMB1_3(:,:,sub) = mean(STFT(:,:,trl1_3(:)),3);
            avMB2_3(:,:,sub) = mean(STFT(:,:,trl2_3(:)),3);
            avMB11_3(:,:,sub) = mean(STFT(:,:,trl11_3(:)),3);% wrong-left
            avMB22_3(:,:,sub) = mean(STFT(:,:,trl22_3(:)),3);% wrong-right

        end

        save(['G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\Motor MuBeta\RT_split3\' fbandname '_' epoch '_CW.mat'], 'avMB1_1', 'avMB2_1', 'avMB11_1', 'avMB22_1', 'avMB1_2', 'avMB2_2', 'avMB11_2', 'avMB22_2', 'avMB1_3', 'avMB2_3', 'avMB11_3', 'avMB22_3');
        clear avMB1_1 avMB2_1 avMB11_1 avMB22_1 avMB1_2 avMB2_2 avMB11_2 avMB22_2 avMB1_3 avMB2_3 avMB11_3 avMB22_3
    end
end
%%%%手动命名保存
% 保存 avMB1 和 avMB2 变量到指定文件路径
%名称规则 波段_切分方法_分类方法_肌肉/条件

%save(['G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\Motor MuBeta\MuBeta_RL_FDIBCP.mat'], 'avMB1', 'avMB2', 'avMB11', 'avMB22');

%           %% Now let's look at Mu/Beta 'MB'
%
% % Let's first verify that there is lateralisation of MB, where the amplitude is lower contralateral to the button that was pressed. Plotting
% % this topography also serves to highlight which electrodes might be best for measuring MB (although note there are some tasks like the
% % continuous dots (2013 paper) and any of our delayed-response tasks, where there is precious little difference in MB amplitude contra/ipsi
% % to responding hand, i.e. not much lateralisation):
%
% % trange = find(Ts>=400,1); % pick a time just before response
% % figure;
% % topoplot(double(avMB1(1:128,trange)-avMB2(1:128,trange)),EEG.chanlocs,'electrodes','labels','colormap','jet');%,'maplimits',scale));
% % %topoplot(double(mean(mean(mean(avMB(1:nchan,trange,1,:,sbj),5),4),2)-mean(mean(mean(avMBr(1:nchan,trange,2,:,sbj),5),4),2)),chanlocs,'electrodes','labels','colormap','jet');%,'maplimits',scale) % some additional arguments that you can add: ,'maplimits',[-1 1]*4 sets the colorscale
% % %topoplot(ERP_right(:, time_indices(i)), EEG_selected_right.chanlocs, 'maplimits', [-max(abs(ERP_right(:))), max(abs(ERP_right(:)))]);
% %
% % title('Left minus Right press')
% % colorbar
%
% num_plots = 15;
% time_points = linspace(-188*6,188*0.75, num_plots);  % in milliseconds
%
% % Find the indices corresponding to the defined time points
% time_indices = arrayfun(@(t) find(Ts >= t, 1), time_points);
%
%
% %左右单独画
% figure;
% for i = 1:num_plots
%     subplot(3, 5, i);  % Create a 2x5 subplot
%     topoplot(double(  mean(avMB1(1:128, time_indices(i),:),3)  -  mean(avMB2(1:128, time_indices(i),:),3)  ) , ...
%          EEG.chanlocs,'electrodes', 'numbers');%,  'maplimits', 0.75*[-0.8,0.8])%,'electrodes', 'labels');
%     % topoplot(double(  avMB2(1:128, time_indices(i)) -  avMB2(1:128, 33) ) , ...
%     %      EEG.chanlocs);%,  'maplimits', 0.75*[-0.8,0.8])%,'electrodes', 'labels');
%      title([num2str(Ts(time_indices(i))), ' ms']);
%     colorbar;
% end
% % %% right only
% % time_indices = find(Ts >= -100 & Ts <= 1200);  % Change the range as required
% % n_subplots = 12;  % Number of subplots
% % time_points = linspace(-100, 1200, n_subplots);  % Select 12 time points between -100 and 1200 ms
% %
% % figure;
% % for i = 1:n_subplots
% %     % Find the closest index in the time vector for the current time point
% %     trange = find(Ts >= time_points(i), 1);
% %
% %     subplot(3, 4, i);  % Create a 3x4 grid of subplots
% %     % Plot the difference between avMB1 and avMB2 at the current time point
% %     topoplot(double(avMB1(1:128, trange)), EEG.chanlocs, ...
% %         'electrodes', 'labels', 'colormap', 'jet');
% %
% %     title(sprintf('Time = %d ms', round(Ts(trange))));  % Title with the current time point
% %     colorbar;
% % end
%
% %% Mu/beta waveforms
%
% ch = [3*32+19 32+22]; % select left/right channels - typically D19 and B22, and that's the case here
%  %ch = [96 76];
% %ch = [108 63];%前 FC34
% %ch = [116 55];%外 C34+
% %ch = [114 53];%里 C34-
% %
% % if sub==6 & epoch == exp.epochs{1}% 为特殊个体选定特殊值
% %
% % ch = [108 54];
% % end
% figure; hold on
%
% plot(Ts,(mean(avMB1(ch(2),:,:),3)+mean(avMB2(ch(1),:,:),3))/2,'r'); % contralateral to correct side
% plot(Ts,(mean(avMB2(ch(2),:,:),3)+mean(avMB1(ch(1),:,:),3))/2,'--b');%(mean(avMB(ch(1),:,1,c,sbj),5)+mean(avMB(ch(2),:,2,c,sbj),5))/2,'--','Color',colours(c,:),'LineWidth',2) % ipsilateral (dashed)
%
% set(gca,'Ydir','reverse') % we often turn the y axis upside down, to show increasing motor preparation (which is reflected in decreasing MB amplitude)
% title('MB: contralateral(red) vs ipsilateral(blue)');
% % % Lateralisation - set it up so upwards means more preparation for the correct alternative
% % figure; hold on
% %
% % plot(Ts,(avMB1(ch(2),:)+avMB1(ch(1),:))/2-(avMB2(ch(1),:)+avMB2(ch(2),:))/2)%,'Color',colours(c,:),'LineWidth',2)
% % title('MB: Lateralisation L-R');
% %
% % % % Response-locked:
% % figure; hold on
% % for c=1:2
% %     plot(Tr,(mean(avMBr(ch(2),:,1,c,sbj),5)+mean(avMBr(ch(1),:,2,c,sbj),5))/2,'Color',colours(c,:),'LineWidth',2) % contralateral to SIDE OF RESPONSE
% %     plot(Tr,(mean(avMBr(ch(1),:,1,c,sbj),5)+mean(avMBr(ch(2),:,2,c,sbj),5))/2,'--','Color',colours(c,:),'LineWidth',2) % ipsilateral (dashed)
% % end
% % set(gca,'Ydir','reverse') % we often turn the y axis upside down, to show increasing motor preparation (which is reflected in decreasing MB amplitude)
% %
% % % notice that the traces coalesce at the time of response, as we've reported lots of times.   