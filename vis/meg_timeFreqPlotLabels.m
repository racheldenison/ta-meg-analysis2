function rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes)

if nargin < 5
    eventTimes = [];
end

set(gca,'YDir','normal')
set(gca,'XTick',xtick)
set(gca,'XTickLabel',toi(xtick))
set(gca,'YTick',ytick)
set(gca,'YTickLabel',foi(ytick))
hold on
for iEv = 1:numel(eventTimes)
    xline(eventTimes(iEv),'k','LineWidth',2);
end