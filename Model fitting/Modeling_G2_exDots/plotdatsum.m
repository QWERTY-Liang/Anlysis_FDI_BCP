function plotdatsum(datsum)
figure; hold on
% colors for 4 coherences:
jt = parula; % a standard colour pallette
col = jt(fliplr(1:21:end),:);
lw = [1 3]; % line widths - error thin and correct thick
% We'll plot it in cumulative form, which is often done in the literature

for ce=1:2 % error and correct
    for c=1:4
        plot((datsum.q(1:end-1,ce,c)+datsum.q(2:end,ce,c))/2,cumsum(datsum.qp(:,ce,c)),'.-','Color',col(c,:),'LineWidth' ,lw(ce))% using the midpoint of the bins as the x-axis
    end
end
cohlevels = [1.07 1.14 2.07 2.14];
legend(num2str(cohlevels'))
end