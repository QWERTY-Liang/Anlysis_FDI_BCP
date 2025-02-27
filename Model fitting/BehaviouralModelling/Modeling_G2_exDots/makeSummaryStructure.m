function datsum = makeSummaryStructure(data,qps,scol)

% gets quantiles (given by qps) and other summary statistics of the raw, single-trial behavioural data in 'data'. 
% Assumes first column is condition, second is outcome (correct/error/miss), third is RT in sec. 
% puts data in a handy structure datsum with these fields:
% n contains the trial numbers for each condition (columns) and outcome, separating error (row 1), correct (row 2) and miss (row 3)
% q contains the quantile RT values (e.g. the value for the .1 quantile is the timepoint below which 10% of the RTs lie)
% qp contains the proportion of trials of a given condition where RT lands in each quantile bin

maxRT = 2.5; % I'm setting this as the upper limit of the slowest quantile bin across the board. If, alternatively, we set this upper limit to be max(data(trl,3)), 
% we would have a gap with no responses between max(data(trl,3)) and 1.5 sec, but yet have some misses beyond that, and
% when simulating with many trials, the model will have no way to reproduce this gap.

% if subject column scol>0, it averages the above 3 fields across subjects indicated in that data column.
% If data not to be averaged, input 0 or [] for scol.
conditions = unique(data(:,1))';
if isempty(scol) | scol==0, % two ways of indicating these data are not to be averaged
    data = [data ones(size(data,1))]; % add a column to say it's all 'subject 1' - this is just so that the same code below can work for either single-subject data or grand averaging across multiple subjects
    scol = size(data,2);
end
subjects = unique(data(:,scol))';

for c=1:length(conditions)
    for s = 1:length(subjects)
        numtr_cond = length(find(data(:,1)==conditions(c) & data(:,scol) == subjects(s))); % number of trials in total for this condition
        for ce = 0:2 % error (0) / correct (1) / miss (2) - for a given condition/subject, these three categories together should make up the whole total of trials
            trl = find(data(:,1)==conditions(c) & data(:,2)==ce & data(:,scol) == subjects(s)); % pull out trials for this condition and subject
            datsum.n(ce+1,c,s) = length(trl); % number of such trials
            if ce<2 % not miss
                % get RT quantiles:
                if datsum.n(ce+1,c,s) % if there are >0 trials with this outcome
                    datsum.q(:,ce+1,c,s) = [min(data(trl,3))-eps prctile(data(trl,3),qps*100) maxRT]'; % The -eps is to shift the minimum RT down slightly just so the fastest RT is included when we use the strict inequality > below
                    % now get proportions in the quantile bins, in exactly the same way as for predicted data in the G^2 function:
                    for qi=1:length(qps)+1  % quantiles
                        datsum.qp(qi,ce+1,c,s) = length(find(data(:,1)==conditions(c) & data(:,2)==ce  & data(:,scol) == subjects(s) & data(:,3)>datsum.q(qi,ce+1,c,s) & data(:,3)<=datsum.q(qi+1,ce+1,c,s)))./numtr_cond;
                    end
                else % if there are no trials (could happen for errors on easy conditions), that subject can't contribute to the group's quantile estimates
                    datsum.q(:,ce+1,c,s) = nan(length(qps)+2,1);
                    datsum.qp(1:length(qps)+1,ce+1,c,s) = 0;
                end
                
            end 
        end
    end
end

% now average across subjects (if this is just one subject it does nothing):
datsum.n = mean(datsum.n, 3);
datsum.q = nanmean(datsum.q, 4);
datsum.qp = mean(datsum.qp, 4);

% Now if this was just one individual, and there were no error trials for a condition, then there would be no quantiles
% to work with for that case. In case that happens, just inherit the same quantiles as for correct trials:
noer = find(isnan(squeeze(datsum.q(1,1,:))));
for c=1:length(noer)
    datsum.q(:,1,noer(c))=datsum.q(:,2,noer(c));
end