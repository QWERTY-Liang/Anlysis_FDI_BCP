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
            % left_channel = 115;  % C3通道号
            % right_channel = 54;
            selected_trials_1 = 999*ones(length(AllBehaviour),3);  % left C4-C3
            selected_trials_2 = 999*ones(length(AllBehaviour),3);  % right C3-C4


            for RTsplit=1:3
                % selected_trials_1 = 999*ones(length(AllBehaviour),RTsplit);  % left C4-C3
                % selected_trials_2 = 999*ones(length(AllBehaviour),RTsplit);  % right C3-C4
                count=0;
                % switch%%%%%%%%%%calculate correct first
                %     case 'CW_SL_b' %correct vs wrong; pre cue SL
                for i=1:length(AllBehaviour)
                    %select trials
                    if AllBehaviour(i,4)==4%AllBehaviour(i,3)==0.14 %AllBehaviour(i,9)==2% &&AllBehaviour(i,4)==4% &&&& AllBehaviour(i,9)==1%&& AllBehaviour(i,3)==0.14%&& AllBehaviour(i,4)==4% && AllBehaviour(i,4)==44
                        selected_trials_1(i,RTsplit)=1;
                        selected_trials_2(i,RTsplit)=0;
                    elseif AllBehaviour(i,4)==44%AllBehaviour(i,3)==0.07%&& AllBehaviour(i,9)==1%&&AllBehaviour(i,4)==4%&& AllBehaviour(i,3)==0.14 %&& AllBehaviour(i,9)==2%&& AllBehaviour(i,3)==0.14%&& AllBehaviour(i,4)==4%&& AllBehaviour(i,4)==44
                        selected_trials_1(i,RTsplit)=0;
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

time_vector = EEG.times;% 保存

selected_idx1 = find(selected_trials_1(:,RTsplit) == 1);
EEG_selected1 = pop_select(EEG, 'trial', selected_idx1);
            erp1_1(:,:,sub,RTsplit) = mean(EEG_selected1.data, 3);% 保存
selected_idx2 = find(selected_trials_2(:,RTsplit) == 1);
EEG_selected2 = pop_select(EEG, 'trial', selected_idx2);
            erp2_1(:,:,sub,RTsplit) = mean(EEG_selected2.data, 3);% 保存
                
            end

           
        end 

  save(['G:\My Drive\Phd\Stage1\BCPvsFDI\E1data_polit\Anlysis_v2\ERP analysis\' 'CPP19_' epoch '_FDIBCP.mat'], 'erp2_1', 'erp1_1', 'time_vector');
       clear erp2_1 erp1_1 time_vector




    end


%% plot ERPs for all RT
% channel_to_plot = 19;
% 
% figure;
% plot(time_vector, mean(mean(erp1_1(channel_to_plot, :,:,:),4),3)   );
% hold on;
% plot(time_vector, mean(mean(erp2_1(channel_to_plot, :,:,:),4),3)  );
% legend('High(0.14)','Low(0.07)')
% xlabel('Time (ms)');
% ylabel('Amplitude (µV)');
% title(['ERP for Channel ', num2str(channel_to_plot)]);
% hold off;
% % Add vertical lines at specified time points
% hold on;
% xline([-1300 -1200, -600, 0, 800, 1500,2000], '--r', {'fix point', 'baseline on', 'cue change', 'evidence on', 'minEvd0.8', 'DDL-1.5s','DDL-2s'});
% hold off;
   