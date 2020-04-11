% meg_plotDecodeGroup.m

%% setup
expt = 'TANoise';
user = 'mcq';
paramType = 'Analysis';

loadAnalysisFiles = 1;

% directories
exptDir = meg_pathToTAMEG(expt, user);
% figDir = sprintf('%s/figures/%s', dataDir, analStr);

% sessions
sessionNames = meg_sessions(expt);
sessionNames = sessionNames([1:8 11:20])

% params
p = meg_params(sprintf('%s_%s', expt, paramType));

%% gather As to make groupA
if loadAnalysisFiles
    groupA = [];
    for iS = 1:numel(sessionNames)
        sessionDir = sessionNames{iS};
        disp(sessionDir)
        
        matDir = sprintf('%s/%s/mat', exptDir, sessionDir);
        dataFile = dir(sprintf('%s/classAcc*nCh50.mat', matDir));
        load(sprintf('%s/%s', matDir, dataFile.name))
        
        if exist('allA','var')
            for iR = 1:length(allA)
                groupA{iS,iR} = allA{iR};
            end
        else
            groupA{iS,1} = A;
        end
    end
end

%% gather data
[nSessions, nRuns] = size(groupA);
nT = numel(groupA{1});
% condNames = {'cueT1','cueT2','neutral'};
condNames = {'cueT1','cueT2'};
% condNames = {'all'};
nCond = numel(condNames);
A = groupA{1,1}(1).(condNames{1});
t = A.classTimes - p.eventTimes(2);

groupData = [];
for iS = 1:nSessions
    for iR = 1:nRuns
        for iT = 1:nT
            for iC = 1:nCond
                condName = condNames{iC};
                groupData(:,iC,iT,iS,iR) = groupA{iS,iR}(iT).(condName).classAcc;
            end
        end
    end
end
        
groupMean = mean(groupData,4);
groupSte = std(groupData,0,4)./sqrt(nSessions);

if strcmp(condNames{1},'cueT1') && strcmp(condNames{2},'cueT2')
    groupDiffData = squeeze(groupData(:,1,:,:) - groupData(:,2,:,:));
    groupDiffData(:,2,:) = -groupDiffData(:,2,:);
    groupDiffMean = mean(groupDiffData,3);
    groupDiffSte = std(groupDiffData,0,3)./sqrt(nSessions);
end

%% plot conds
figure
colors = get(gca,'ColorOrder');
close(gcf)
% colors = p.cueColors;

ylims = [45 65];
figure('Position',[200 500 900 400])
p1 = [];
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

%% plot runs
groupData = squeeze(groupData);
groupMean = squeeze(groupMean);
groupSte = squeeze(groupSte);
runNames = num2str(nTopCh');

figure
colors = get(gca,'ColorOrder');
close(gcf)
% colors = p.cueColors;

ylims = [45 70];
figure('Position',[200 500 900 400])
for iT = 1:nT
    subplot(1,nT,iT)
    hold on
    for iR = 1:nRuns
%         shadedErrorBar(t, groupMean(:,iT,iR), groupSte(:,iT,iR), {'Color', colors(iR,:)}, 1)
        p1(iR) = plot(t, groupMean(:,iT,iR), 'Color', colors(iR,:));
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
legend(p1, runNames)

% individual data
ylims = [40 max(groupData(:))*1.1];
for iT = 1:nT
    figure
    for iS = 1:nSessions
        subplot(10,2,iS)
        hold on
        plot(t, squeeze(groupData(:,iT,iS,:)));
        plot(t([1 end]),[50 50],'k')
        xlim(t([1 end]))
        ylim(ylims)
        vline(0,'--k')
        if iS==1
            xlabel('Time (ms)')
            ylabel('Classification accuracy (%)')
            legend(runNames)
        end
    end
    rd_supertitle2(sprintf('T%d',iT))
end

% scatter plots
groupDataVals.mean = squeeze(mean(groupData,1));
groupDataVals.max = squeeze(max(groupData,[],1));

lims.mean = [40 max(groupDataVals.mean(:))*1.1];
lims.max = [40 max(groupDataVals.max(:))*1.1];

measureNames = fields(groupDataVals);

for iT = 1:nT
    figure('Position',[800 500 900 400])
    for iM = 1:numel(measureNames)
        m = measureNames{iM};
        for iR = 1:nRuns-1
            subplot(2,nRuns-1,(nRuns-1)*(iM-1)+iR)
            hold on
            plot(lims.(m), lims.(m), 'k')
            plot(squeeze(groupDataVals.(m)(iT,:,1)), squeeze(groupDataVals.(m)(iT,:,iR+1)),...
                '.','Color',colors(iR+1,:),'MarkerSize', 30)
            xlim(lims.(m))
            ylim(lims.(m))
            axis square
            title(runNames(iR+1,:))
            if iR==1 && iM==2
                xlabel('all channels')
                ylabel('subset of channels')
            end
            if iR == nRuns-1
                yyaxis right
                
                set(gca,'YTick',[],'YColor','w')
                ylabel(m,'Color','k')
            end
        end
    end
    rd_supertitle2(sprintf('T%d classification accuracy, all trials',iT))
end

runIdx = 3;
figure
for iM = 1:numel(measureNames)
    m = measureNames{iM};
    vals = groupDataVals.(m)(:,:,runIdx) - groupDataVals.(m)(:,:,1);
    
    subplot(1,numel(measureNames),iM)
    hold on
    plot(vals')
    plot([0 20], [0 0], 'k')
    legend(targetNames)
    xlabel('session')
    ylabel('change in decoding from 157 ch')
    title(sprintf('%s, %s ch', m, runNames(runIdx,:)))
end

meanAve = squeeze(mean(mean(groupDataVals.mean,1),2))';
maxAve = squeeze(mean(mean(groupDataVals.max,1),2))';
disp([nTopCh; meanAve; maxAve])
