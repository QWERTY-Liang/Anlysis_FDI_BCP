function plotdatsum(datsum)
figure; hold on
% colors for 4 coherences:
% jt = parula; % a standard colour pallette
% col = jt(fliplr(1:21:end),:);
col = lines(4);
lw = [1 3]; % line widths - error thin and correct thick
% We'll plot it in cumulative form, which is often done in the literature

for ce=1:2 % error and correct
    for c=1:4
        plot((datsum.q(1:end-1,ce,c)+datsum.q(2:end,ce,c))/2,cumsum(datsum.qp(:,ce,c)),'.-','Color',col(c,:),'LineWidth' ,lw(ce))% using the midpoint of the bins as the x-axis
    end
end
% cohlevels = [6 14 26 48];
% legend(num2str(cohlevels'))
cohlevels = ["FDI_low", "FDI_high", "BCP_low", "BCP_high"];
legend(cohlevels);
end