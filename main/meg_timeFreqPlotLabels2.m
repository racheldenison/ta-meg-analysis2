function meg_timeFreqPlotLabels2(toi,foi,xtick,ytick,eventTimes)

if nargin < 5
    eventTimes = [];
end

set(gca,'YDir','normal')
set(gca,'XTick',xtick)
set(gca,'XTickLabel',xtick)
set(gca,'YTick',ytick)
set(gca,'YTickLabel',ytick)
hold on
for iEv = 1:numel(eventTimes)
    xline(eventTimes(iEv),'k','LineWidth',2);
end