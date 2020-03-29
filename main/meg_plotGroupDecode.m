% meg_plotGroupDecode.m

p = meg_params('TA2_Analysis');

nSessions = numel(groupA);
nT = 2;
% condNames = {'cueT1','cueT2','neutral'};
condNames = {'cueT1','cueT2'};
nCond = numel(condNames);
A = groupA{1}(1).cueT1;
t = A.classTimes - p.eventTimes(2);

groupData = [];
for iS = 1:nSessions
    for iT = 1:nT
        for iC = 1:nCond
            condName = condNames{iC};
            groupData(:,iC,iT,iS) = groupA{iS}(iT).(condName).classAcc;
        end
    end
end
        
groupMean = mean(groupData,4);
groupSte = std(groupData,0,4)./sqrt(nSessions);

groupDiffData = squeeze(groupData(:,1,:,:) - groupData(:,2,:,:));
groupDiffData(:,2,:) = -groupDiffData(:,2,:);
groupDiffMean = mean(groupDiffData,3);
groupDiffSte = std(groupDiffData,0,3)./sqrt(nSessions);

%% plot
figure
colors = get(gca,'ColorOrder');
close(gcf)
% colors = p.cueColors;

ylims = [45 65];
figure('Position',[200 500 900 400])
for iT = 1:nT
    subplot(1,nT,iT)
    hold on
    for iC = 1:nCond
        shadedErrorBar(t, groupMean(:,iC,iT), groupSte(:,iC,iT), {'Color', colors(iC,:)}, 1)
%         shadedErrorBar(t, groupMean(:,iC,iT), groupDiffSte(:,iT), {'Color', colors(iC,:)}, 1)
        p1(iC) = plot(t, groupMean(:,iC,iT), 'Color', colors(iC,:));
    end
    plot(t([1 end]),[50 50],'k')
    title(sprintf('T%d',iT))
    xlim(t([1 end]))
    ylim(ylims)
    vline(0,'--k')
    if iT==1
        xlabel('Time (ms)')
        ylabel('Classification accuracy (%)')
    end
end
legend(p1, condNames)

ylims = [-10 15];
figure
for iT = 1:nT
    subplot(nT,1,iT)
    hold on
    shadedErrorBar(t, groupDiffMean(:,iT), groupDiffSte(:,iT), 'k', 1)
    plot(t([1 end]),[0 0],'k')
    title(sprintf('T%d',iT))
    ylim(ylims)
    if iT==nT
        xlabel('Time (ms)')
        ylabel('Classification accuracy, cued - uncued (%)')
    end
end