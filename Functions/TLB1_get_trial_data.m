function [Contrast Muscle Perf Rt EVtime TotEMG TotEMG_BCP RespLR CorrLR FirstTilt] = TLB1_get_trial_data(filename,exp)


load([exp.behpath  filename]);


Contrast = par.deltaC(par.cond);%0.07/0.14
Muscle=par.muscleorder;% 1=FDI 2-BCP
Perf = perf;%(1=correct, 2=error, 4=no response, 3 too early, 6 slow, 5 wrong muscle)
Rt = ThresholdTime - EvOn_time;
EVtime= evoff - EvOn_time;% 
TotEMG = [trialBaseline; trialEMG; data_postEMG];
TotEMG_BCP = [trialBaseline_becips; trialEMG_becips; data_postEMG_becips];
RespLR = respLR;
CorrLR = par.LR;
FirstTilt = par.firstTilt;
% Duration = par.evdur(par.cond);
% Condition = par.cond;
% %Threshold = trial_thresh;
% LenBase = length(trialBaseline);
% LenTrial = length(trialEMG);
% LenPost = length(data_postEMG);


end
