% meg_plotGroupDecode.m

nSessions = numel(groupA);
nT = 2;
% condNames = {'cueT1','cueT2','neutral'};
condNames = {'cueT1','cueT2'};
nCond = numel(condNames);
A = groupA{1}(1).cueT1;
t = A.classTimes - 1000;

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

groupDiffData = squeeze(groupData(:,1,:,:) - groupData(:,2,:,:));
groupDiffData(:,2,:) = -groupDiffData(:,2,:);
groupDiffMean = mean(groupDiffData,3);
groupDiffSte = std(groupDiffData,0,3)./sqrt(nSessions);

%% plot
ylims = [45 65];
figure
for iT = 1:nT
    subplot(1,nT,iT)
    hold on
    plot(t, groupMean(:,:,iT))
    plot(t([1 end]),[50 50],'k')
    title(sprintf('T%d',iT))
    ylim(ylims)
    if iT==1
        xlabel('Time (ms)')
        ylabel('Classification accuracy (%)')
    end
end
legend(condNames)

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