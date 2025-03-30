function g2 = G2(pvar,pfix,Sel,datsum,qps,modelfn,N,seed)

% G^2 is a goodness of fit metric we're using as the objective function to be minimised, which is calculated as the sum
% over a whole bunch of no_c*po*log(po/pp) terms, where no_c is the total number of observed trials in a condition c, po is the
% observed proportion of trials in that condition with a certain outcome category defined as a quantile bin in either
% the correct or error disribution, or misses, and pp is the equivalent in the predicted (model-simulated) data

eval(['simdata = ' modelfn '(pvar,pfix, Sel,N,seed);'])

% Get the predicted proportions pp of simulated trials falling in each of the outcome categories defined for the real data - the quantile bins for corrects, errors, and the misses. 
pp=[]; po=[]; % predicted and observed proportions
no_c=[]; % total number of trials in the OBSERVED data in the corresponding condition. It's in datsum but we pull it out in the loop below to gather the correct n values alongside the predicted proportions for the final g2 calculation at the end
for c=1:size(datsum.n,2)   %  conditions
    np_c=length(find(simdata(:,1)==c));   % total number of trials in the PREDICTED data for the current condition
    pp_c=[]; % predicted proportions for just this condition
    for ce=0:1 % correct/error
        for qi=1:length(qps)+1  % quantiles
            pp_c = [pp_c length(find(simdata(:,1)==c & simdata(:,2)==ce & simdata(:,3)>datsum.q(qi,ce+1,c) & simdata(:,3)<=datsum.q(qi+1,ce+1,c)))./np_c];
            po = [po datsum.qp(qi,ce+1,c)]; % corresponding observed proportions in the real data
            no_c = [no_c sum(datsum.n(:,c))];
        end
    end
    % Misses:
    pp_c = [pp_c length(find(simdata(:,1)==c & simdata(:,2)==2))./np_c]; % the simulation function sets to 2 the outcome of any trial where no bound is crossed anytime up to maxRT of 1.5 sec
    po = [po 1-sum(sum(datsum.qp(:,:,c)))]; % corresponding observed proportions of misses in the real data. In THAT case, because we know anything not captured as an RT must be a miss, we can compute it as the proportion left over. This is not true of predicted pp because RTs could be generated outside the observed data range within which we compare the proportions. See below
    no_c = [no_c sum(datsum.n(:,c))]; % number of trials in this condition. This will be the same number repeated for all outcomes within that condition - this is for convenience of the elementwise operations done on the vectors in the last line of code
    
    % In the n_c*po*log(po/pp) terms, we'll get values of Inf(inity) if pp=0, which may hamper the algorithm distinguishing where to go, so
    % we apply a correction to avoid it, and re-normalise across all predicted data outcomes so the pp's (including the residual proportion
    % landing completely outside the range of observed data, 'pp_out') still sum to 1
    pp_out = 1-sum(pp_c);
    pp_c(pp_c<0.0001) = 0.0001; % only apply correction where needed
    pp_c = pp_c./(sum(pp_c)+pp_out); % re-normalise so all predicted proportions sum to 1 again, including outside values. Won't change anything if there are no pp_c<0.0001.
    % Note that pp_c is used in the below g2 calc whereas pp_out trials are just thrown away. The proportions in pp_c by itself won't sum to 1 in the event p_out>0, 
    % but this is as it should be - the denominator for all proportions must always be the total trials of that condition that were
    % simulated, and we have to allow for the fact that some RTs in the predicted data may land in the 'dead zone' RT<minRealRT where there are no real RTs

    % now append to overall pp vector:
    pp = [pp pp_c];
end

% compute G-square:
g2=2*nansum(no_c.*po.*log(po./pp));
% I use nansum here because if po = 0 for any outcome of the observed data, 0*log(0/x) will be NaN which presumably won't help the search
% algorithm. In this particular example, where we have >0 trials in every outcome cell in the observed data we have no case of po=0, but you
% might have this in other experiments