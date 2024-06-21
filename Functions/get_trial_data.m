function [Contrast Duration Condition Threshold RespLR CorrLR Rt LenBase LenTrial LenPost TotEMG Perf] = get_trial_data(filename)
load(filename);
Contrast = par.deltaC(par.cond);
Duration = par.evdur(par.cond);
Condition = par.cond;
Threshold = trial_thresh;
RespLR = respLR;
CorrLR = par.LR;
Rt = ThresholdTime - EvOn_time;
LenBase = length(trialBaseline);
LenTrial = length(trialEMG);
LenPost = length(data_postEMG);
TotEMG = [trialBaseline; trialEMG; data_postEMG];
Perf = perf;
%FirstTilt = par.firstTilt;
end
