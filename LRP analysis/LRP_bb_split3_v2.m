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

%For left-hand responses, C4 is the contralateral (right hemisphere) signal and C3 is the ipsilateral (left hemisphere) signal.
%For right-hand responses, C3 is the contralateral (left hemisphere) signal and C4 is the ipsilateral (right hemisphere) signal.
%LRP=Contra−Ipsi
%%条件1
            % left_channel = 115;  % C3通道号
            % right_channel = 54;
            selected_trials_1 = 999*ones(length(AllBehaviour),3);  % left C4-C3
            selected_trials_2 = 999*ones(length(AllBehaviour),3);  % right C3-C4

%%条件2
            % left_channel = 115;  % C3通道号
            % right_channel = 54;
            selected_trials2_1 = 999*ones(length(AllBehaviour),3);  % left C4-C3
            selected_trials2_2 = 999*ones(length(AllBehaviour),3);  % right C3-C4 


            for RTsplit=1:3
                % selected_trials_1 = 999*ones(length(AllBehaviour),RTsplit);  % left C4-C3
                % selected_trials_2 = 999*ones(length(AllBehaviour),RTsplit);  % right C3-C4
                count=0;
                % switch%%%%%%%%%%calculate correct first
                %     case 'CW_SL_b' %correct vs wrong; pre cue SL
                for i=1:length(AllBehaviour)
 %条件一
                    if AllBehaviour(i,8)==1&&AllBehaviour(i,4)==4% && AllBehaviour(i,3)==0.14%AllBehaviour(i,3)==0.14 %AllBehaviour(i,9)==2% &&AllBehaviour(i,4)==4% &&&& AllBehaviour(i,9)==1%&& AllBehaviour(i,3)==0.14%&& AllBehaviour(i,4)==4% && AllBehaviour(i,4)==44
                        selected_trials_1(i,RTsplit)=1;
                        selected_trials_2(i,RTsplit)=0;
                    elseif AllBehaviour(i,8)==2 &&AllBehaviour(i,4)==4%&& AllBehaviour(i,3)==0.14%AllBehaviour(i,3)==0.07%&& AllBehaviour(i,9)==1%&&AllBehaviour(i,4)==4%&& AllBehaviour(i,3)==0.14 %&& AllBehaviour(i,9)==2%&& AllBehaviour(i,3)==0.14%&& AllBehaviour(i,4)==4%&& AllBehaviour(i,4)==44
                        selected_trials_1(i,RTsplit)=0;
                        selected_trials_2(i,RTsplit)=1;
                    end
          


                    %exclude invalide(too early or wrong muscle)
                    if AllBehaviour(i,5)==3 || AllBehaviour(i,5)==5|| AllBehaviour(i,5)==4%|| AllBehaviour(i,5)==6
                        selected_trials_1(i,RTsplit)=0;%change back to unselected
                        selected_trials_2(i,RTsplit)=0;
                        count=count+1; %计数丢掉的trial
                    end
%条件2
 if AllBehaviour(i,8)==1 &&AllBehaviour(i,4)==44%&& AllBehaviour(i,3)==0.07%AllBehaviour(i,3)==0.14 %AllBehaviour(i,9)==2% &&AllBehaviour(i,4)==4% &&&& AllBehaviour(i,9)==1%&& AllBehaviour(i,3)==0.14%&& AllBehaviour(i,4)==4% && AllBehaviour(i,4)==44
                        selected_trials2_1(i,RTsplit)=1;
                        selected_trials2_2(i,RTsplit)=0;
                    elseif AllBehaviour(i,8)==2 &&AllBehaviour(i,4)==44%&& AllBehaviour(i,3)==0.07%AllBehaviour(i,3)==0.07%&& AllBehaviour(i,9)==1%&&AllBehaviour(i,4)==4%&& AllBehaviour(i,3)==0.14 %&& AllBehaviour(i,9)==2%&& AllBehaviour(i,3)==0.14%&& AllBehaviour(i,4)==4%&& AllBehaviour(i,4)==44
                        selected_trials2_1(i,RTsplit)=0;
                        selected_trials2_2(i,RTsplit)=1;
                    end
          


                    %exclude invalide(too early or wrong muscle)
                    if AllBehaviour(i,5)==3 || AllBehaviour(i,5)==5|| AllBehaviour(i,5)==4%|| AllBehaviour(i,5)==6
                        selected_trials2_1(i,RTsplit)=0;%change back to unselected
                        selected_trials2_2(i,RTsplit)=0;
                        count=count+1; %计数丢掉的trial
                    end





                    %把不属于时间的剔除
                    if RTsplit==1
                        if RT(i)>=RT33%只保留偏小三等分点下的
                            selected_trials_1(i,RTsplit)=0;%change back to unselected
                            selected_trials_2(i,RTsplit)=0;

                             selected_trials2_1(i,RTsplit)=0;%change back to unselected
                            selected_trials2_2(i,RTsplit)=0;
                        end
                        trl1_1=find(selected_trials_1(:,RTsplit)==1);
                        trl2_1=find(selected_trials_2(:,RTsplit)==1);

                        trl11_1=find(selected_trials2_1(:,RTsplit)==1);
                        trl22_1=find(selected_trials2_2(:,RTsplit)==1);

                    elseif RTsplit==2
                        if RT(i)<=RT33 || RT(i)>=RT66%只保留两个三等分点中的
                            selected_trials_1(i,RTsplit)=0;%change back to unselected
                            selected_trials_2(i,RTsplit)=0;

                            selected_trials2_1(i,RTsplit)=0;%change back to unselected
                            selected_trials2_2(i,RTsplit)=0;
                        end

                        trl1_2=find(selected_trials_1(:,RTsplit)==1);
                        trl2_2=find(selected_trials_2(:,RTsplit)==1);

                        trl11_2=find(selected_trials2_1(:,RTsplit)==1);
                        trl22_2=find(selected_trials2_2(:,RTsplit)==1);
                    elseif RTsplit==3
                        if RT(i)<=RT66%只保留偏da三等分点上的
                            selected_trials_1(i,RTsplit)=0;%change back to unselected
                            selected_trials_2(i,RTsplit)=0;
                           
                             selected_trials2_1(i,RTsplit)=0;%change back to unselected
                            selected_trials2_2(i,RTsplit)=0;
                        end
                        trl1_3=find(selected_trials_1(:,RTsplit)==1);
                        trl2_3=find(selected_trials_2(:,RTsplit)==1);

                        trl11_3=find(selected_trials2_1(:,RTsplit)==1);
                        trl22_3=find(selected_trials2_2(:,RTsplit)==1);

                    end

                end  
            % end

time_vector = EEG.times;% 保存

right_channel=[63,64, 53,54];% 保存

left_channel=[108,109, 114,115];% 保存

%条件1
% Extract selected trials---for correct
selected_idx_1 = find(selected_trials_1(:,RTsplit) == 1);
EEG_selected_left = pop_select(EEG, 'trial', selected_idx_1);
% Compute the ERP
ERP_left(:,:,sub,RTsplit) = mean(EEG_selected_left.data, 3);  % averaging across the third dimension (trials)

LRP_left(:,:,sub,RTsplit) = (mean(ERP_left(right_channel, :,sub,RTsplit),1) - mean(ERP_left(left_channel, :,sub,RTsplit),1)) / 2;%这里先平均再相减

%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
%time_vector_left = EEG_selected_left.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% Extract selected trials---for wrong
selected_idx_right = find(selected_trials_2(:,RTsplit) == 1);
EEG_selected_right = pop_select(EEG, 'trial', selected_idx_right);
% Compute the ERP
ERP_right(:,:,sub,RTsplit) = mean(EEG_selected_right.data, 3);  % averaging across the third dimension (trials)
LRP_right(:,:,sub,RTsplit) = (mean(ERP_right(left_channel, :,sub,RTsplit),1) - mean(ERP_right(right_channel, :,sub,RTsplit),1)) / 2;%这里先平均再相减
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
%time_vector_right = EEG_selected_right.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

LRP_con1(:,:,sub,RTsplit)=LRP_left(:,:,sub,RTsplit)+LRP_right(:,:,sub,RTsplit); %保存

%条件2

% Extract selected trials---for correct
selected_idx2_1 = find(selected_trials2_1(:,RTsplit) == 1);
EEG_selected_left2 = pop_select(EEG, 'trial', selected_idx2_1);
% Compute the ERP
ERP_left2(:,:,sub,RTsplit) = mean(EEG_selected_left2.data, 3);  % averaging across the third dimension (trials)

LRP_left2(:,:,sub,RTsplit) = (mean(ERP_left2(right_channel, :,sub,RTsplit),1) - mean(ERP_left2(left_channel, :,sub,RTsplit),1)) / 2;%这里先平均再相减

%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
%time_vector_left = EEG_selected_left.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

% Extract selected trials---for wrong
selected_idx_right2 = find(selected_trials2_2(:,RTsplit) == 1);
EEG_selected_right2 = pop_select(EEG, 'trial', selected_idx_right2);
% Compute the ERP
ERP_right2(:,:,sub,RTsplit) = mean(EEG_selected_right2.data, 3);  % averaging across the third dimension (trials)
LRP_right2(:,:,sub,RTsplit) = (mean(ERP_right2(left_channel, :,sub,RTsplit),1) - mean(ERP_right2(right_channel, :,sub,RTsplit),1)) / 2;%这里先平均再相减
%ERP = mean(EEG_selected.data, 3);  % averaging across the third dimension (trials)
% Plot the ERP for a specific channel
%time_vector_right = EEG_selected_right.times;  % time points from the EEG structure
% channel_to_plot = 19;  % channel number to plot

LRP_con2(:,:,sub,RTsplit)=LRP_left2(:,:,sub,RTsplit)+LRP_right2(:,:,sub,RTsplit); %保存









            end

           
        end 

   save(['G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\LRP analysis\' 'LRP_' epoch '_FDIBCP.mat'], 'LRP_con1', 'LRP_con2', 'time_vector','right_channel','left_channel','ERP_left2','ERP_left','ERP_right2','ERP_right');
        clear LRP_con1 LRP_con2 time_vector right_channel left_channel ERP_left2 ERP_left ERP_right2 ERP_right




    end
