% April 2020

% setup 
expt = 'TA2'; 
user = 'karen';
slice = 'cue_subsample50';  % cueAcc_correct'; % cueAcc 'cueAcc_incorrect'
exptDir = meg_pathToTAMEG(expt, user);
% analStr = 'bietfp'; 

saveFigs = 1; 

%% setup
% directory
figDir = sprintf('%s/Group/figures/%s_max/',exptDir,slice);
if ~exist(figDir,'dir')
    mkdir(figDir)
end
[sessionNames,subjectNames] = meg_sessions(expt); 

%% setup data and ch matrix w new nomenclature 
for iSession = 1:numel(sessionNames)
    groupC(iSession).channelsRanked = groupC(iSession).sortChByProm;
    groupC(iSession).channelDirection = groupC(iSession).chPromDir;
end

%% extract data of interest 
D = []; 

% cue 
D.cueT1 = groupDavg.cueT1;
D.cueT2 = groupDavg.cueT2;
D.cueN = groupDavg.cueN;

% incorrect 
% D.cueT1 = groupDavg.cueT1_incorrect; 
% D.cueT2 = groupDavg.cueT2_incorrect; 
% D.cueN = groupDavg.cueN_incorrect; 

% correct
% D.cueT1 = groupDavg.cueT1_correct; 
% D.cueT2 = groupDavg.cueT2_correct; 
% D.cueN = groupDavg.cueN_correct; 

%% flip data
D_bietfp_flipped = []; 
for iSession = 1:numel(sessionNames)
    chDir = groupC(iSession).channelDirection; 
    D_bietfp_flipped.cueT1(:,:,iSession) = D.cueT1(:,:,iSession).*chDir'; 
    D_bietfp_flipped.cueT2(:,:,iSession) = D.cueT2(:,:,iSession).*chDir'; 
    if strcmp(expt,'TA2')
        D_bietfp_flipped.cueN(:,:,iSession) = D.cueN(:,:,iSession).*chDir'; 
    end
end

%% subsample checks
if subSample 
    %%  get trial data
    trialD.cueT1 = groupD2.cueT1All;
    trialD.cueT2 = groupD2.cueT2All;
    trialD.cueN = groupD2.neutralAll;
    
    %% 
    nSample = 50;
    fields = fieldnames(groupD{1});
    Dsubsample = []; 
    for iF = 1:numel(fields)
        for i = 1:numel(sessionNames)
        % Dsubsample.nSample = nSample;
        fieldVals = groupD{i}.(fields{iF});
        nTrials = size(fieldVals,3); 
        
        idx = []; subsampleVals = []; avgSubsampleVals = []; avgAvgSubSample = []; 
        for iS = 1:100 % sample 100 times
            idx = randperm(nTrials,nSample);
            subsampleVals = fieldVals(:,:,idx); % time x ch x subsampled trials
            avgSubsampleVals(:,:,iS) = squeeze(nanmean(subsampleVals,3)); %  time x ch x 1 avg trial 
        end 
        avgAvgSubSample = squeeze(nanmean(avgSubsampleVals,3)); 
        Dsubsample.(fields{iF})(:,:,i) = avgAvgSubSample; 
        end
    end
    
    Dsubsample50 = Dsubsample; 
    
    %% flip 
    D = Dsubsample50; 
    D_bietfp_flipped = [];
    for iSession = 1:numel(sessionNames)
        channels = groupC(iSession).channelsRanked; 
        chDir = groupC(iSession).channelDirection(channels);
        D_bietfp_flipped.cueT1(:,:,iSession) = D.cueT1(:,:,iSession).*chDir';
        D_bietfp_flipped.cueT2(:,:,iSession) = D.cueT2(:,:,iSession).*chDir';
        if strcmp(expt,'TA2')
            D_bietfp_flipped.cueN(:,:,iSession) = D.neutral(:,:,iSession).*chDir';
        end
    end
end

%% 
channelGroups = [1,2,3,5,7,10,15,20]; 
% for nChOI = 1:20
for i = 1:numel(channelGroups)
    % nChOI = 1:channelGroups(i); 
    nChOI = 1:5; 
    disp(nChOI)
    [D,A] = meg_peakAnalysis(D_bietfp_flipped,nChOI,groupC,expt,'karen',slice); % meg_plotERFMax_chOI
    close all 
end
% end
    
%% ttest max and max idx by cue
[hT1max, pvalT1max] = ttest(A.T1.cueT1.max,A.T1.cueT2.max,'alpha',0.05);
[hT2max, pvalT2max] = ttest(A.T2.cueT1.max,A.T2.cueT2.max,'alpha',0.05);

[hT1maxIdx, pvalT1maxIdx] = ttest(A.T1.cueT1.maxIdx,A.T1.cueT2.maxIdx,'alpha',0.05);
[hT2maxIdx, pvalT2maxIdx] = ttest(A.T2.cueT1.maxIdx,A.T2.cueT2.maxIdx,'alpha',0.05);

%% groupA 

groupA = [];
for i = [1,2,3,5,7,10,15,20]
    % cond = 'cue_max'; % cueAcc_max
    % groupA.cond = cond; 
        
    switch expt 
        case 'TANoise' 
            groupA.expt = 'TANoise'; 
            groupA.slice = slice; 
            load(sprintf('/Users/kantian/Dropbox/Data/TANoise/MEG/Group/figures/%s/nCh_1_%d/max.mat',cond,i))
        case 'TA2' 
            groupA.expt = 'TA2'; 
            groupA.slice = slice; 
            load(sprintf('%s/nCh_1_%d/max.mat',figDir,i))
    end
    
    newField = sprintf('chTop%d',i);
    vals = A;
    groupA.(newField) = vals;
    
    % paired ttest 
    [hT1max, pvalT1max, ciT1max, statsT1max] = ttest(A.T1.cueT1.max,A.T1.cueT2.max,'alpha',0.05);
    [hT2max, pvalT2max, ciT2max, statsT2max] = ttest(A.T2.cueT1.max,A.T2.cueT2.max,'alpha',0.05);
    
    [hT1maxIdx, pvalT1maxIdx, ciT1maxIdx, statsT1maxIdx] = ttest(A.T1.cueT1.maxIdx,A.T1.cueT2.maxIdx,'alpha',0.05);
    [hT2maxIdx, pvalT2maxIdx, ciT2maxIdx, statsT2maxIdx] = ttest(A.T2.cueT1.maxIdx,A.T2.cueT2.maxIdx,'alpha',0.05);
    
    % save ttest 
    groupA.(newField).hT1max = hT1max; 
    groupA.(newField).pvalT1max = pvalT1max;
    
    groupA.(newField).hT2max = hT2max; 
    groupA.(newField).pvalT2max = pvalT2max;
    
    groupA.(newField).hT1maxIdx = hT1maxIdx; 
    groupA.(newField).pvalT1maxIdx = pvalT1maxIdx;
    
    groupA.(newField).hT2maxIdx = hT2maxIdx; 
    groupA.(newField).pvalT2maxIdx = pvalT2maxIdx;
end

%%
save(sprintf('%s/groupA.mat',figDir),'groupA')

%% setup for figures
p = meg_params(sprintf('%s_Analysis',expt));
tIdx = 1:p.trialTime; % epoch time, 1 indexed
t = p.tstart:p.tstop; % epoch time, relative to precue
targetNames = p.eventNames(2:3); 
[sessionNames,subjectNames] = meg_sessions(expt); 
switch expt 
    case 'TANoise'
        windowPostTarget = 80; % 70 min window post target to start looking for evoked max, 0 if right after target onset 
        windowSize = 170;  % 190 TANoise(1)  window length to look for target evoked max 
    case 'TA2'
        windowPostTarget = 100;
        windowSize = 200;
end
fields = fieldnames(D.subject); 

%% plot channel grouping max
figure
set(gcf,'Position',[100 100 1000 400])

for iT = 1:numel(targetNames)
    target = targetNames(iT);
    idxTargetName = find(strcmp(p.eventNames,target));
    idxT = find(t==p.eventTimes(idxTargetName));
    window = idxT:idxT+windowPostTarget+windowSize;
    for iC = 1:numel(channelGroups)
        channelGroup = sprintf('chTop%d',channelGroups(iC)); 
        for iF = 1:numel(fields)
        max = groupA.(channelGroup).(targetNames{iT}).(fields{iF}).max; 
        maxMean = nanmean(max); 
        subplot (1,numel(targetNames),iT)
        hold on
            errV = errorbar(iC,...%groupA.(channelGroup).(targetNames{iT}).(fields{iF}).maxIdxMean,... % X: max idx
                maxMean,... % Y: max
                nanstd(max)/sqrt(numel(subjectNames)),... % vertical error bar: max ste
                'Color',p.cueColors(iF,:),'LineWidth',2,'CapSize',10);
        end
        for iF = 1:numel(fields)
            p1 = plot(iC,... %groupA.(channelGroup).(targetNames{iT}).(fields{iF}).maxIdxMean,...
                maxMean,...
                'o','MarkerSize',8);
            p1.Color = p.cueColors(iF,:);
            p1.MarkerFaceColor = p.cueColors(iF,:);
        end
        % ylim([0,0.8e-13])
        ylim([0,2.8e-13])
        xlim([0,numel(channelGroups)+1])
        xlabel('Number of channels')
        xticks(1:numel(channelGroups))
        xticklabels(string(num2cell(channelGroups)))
        ylabel('Amplitude (T)')
        vline(p.eventTimes(iT+1),'k:','target')
        title(sprintf('%s max by cue condition',targetNames{iT}))
    end
    % disp(pval(iC)) 
end
if saveFigs 
    saveas(gcf,sprintf('%s/maxAmp.svg',figDir)) 
end

%% plot channel grouping selection max idx

figure
set(gcf,'Position',[100 100 1000 200])

for iT = 1:numel(targetNames)
    target = targetNames(iT);
    idxTargetName = find(strcmp(p.eventNames,target));
    idxT = find(t==p.eventTimes(idxTargetName));
    window = idxT:idxT+windowPostTarget+windowSize;
    for iC = 1:numel(channelGroups)
        channelGroup = sprintf('chTop%d',channelGroups(iC)); 
        subplot (1,numel(targetNames),iT)
        hold on
        for iF = 1:numel(fields)
            maxIdx = groupA.(channelGroup).(targetNames{iT}).(fields{iF}).maxIdx; 
            maxMean = nanmean(maxIdx); 
            subplot (1,numel(targetNames),iT)
            % vertical error bar
            errV = errorbar(iC,...
                maxMean-p.eventTimes(iT+1),... 
                nanstd(maxIdx)/sqrt(numel(subjectNames)),... % vertical error bar: max ste
                'Color',p.cueColors(iF,:),'LineWidth',2,'CapSize',10);
            p1 = plot(iC,...
                maxMean-p.eventTimes(iT+1),...
                'o','MarkerSize',8);
            p1.Color = p.cueColors(iF,:);
            p1.MarkerFaceColor = p.cueColors(iF,:);
        end
        ylim([90,150])
        xlim([0,numel(channelGroups)+1])
        xticks(1:numel(channelGroups))
        xticklabels(string(num2cell(channelGroups)))
        xlabel('Number of channels')
        ylabel('Time post target (ms)')
        title(sprintf('%s Peak latency',targetNames{iT}))
    end
end
if saveFigs 
    saveas(gcf,sprintf('%s/maxIdx.svg',figDir)) 
end

%% if sig, plot to check cue max by subject
figure
hold on 
plot(A.T1.cueT1.maxIdx)
plot(A.T2.cueT1.maxIdx)
title('ch 1:2')

%% setup new max abs amp stuff 
% for iSession = 1:numel(sessionNames)
%     chDir = []; 
%     chDir = groupC_maxAbs(iSession).chDir; 
%     D_input.cueT1(:,:,iSession) = D_bietfp.cueT1(:,:,iSession).*chDir;
%     D_input.cueT2(:,:,iSession) = D_bietfp.cueT2(:,:,iSession).*chDir;
%     D_input.cueN(:,:,iSession)  = D_bietfp.cueN(:,:,iSession).*chDir;
% end

%% get rid of wonky 1 for R1547 session 2, session20 

% groupC2 = groupC;
% % delete 1st channel 
% % groupC2(20).channelsRanked(1) = []; 
% % groupC2(20).dirRanked(1) = []; 
% groupC2(9).dirRanked(1) = groupC(9).dirRanked(1)*(-1); 
% groupC2(14).dirRanked(1) = groupC(14).dirRanked(1)*(-1); 
% 
% groupC2(9).dirRanked(1) = groupC(9).dirRanked(1)*(-1); 
% groupC2(14).dirRanked(1) = groupC(14).dirRanked(1)*(-1); 
% 
% groupC2(9).chPromDir(1) = groupC2(9).dirRanked(1) 
% groupC2(14).chPromDir(111) = groupC2(14).dirRanked(1); 
% 
% %% unflip data, so raw 
% for iSession = 1:numel(sessionNames)
%     for iC = 1:157
%         if groupC(iSession).chPromDir(iC) == -1
%             D6.cueT1(:,iC,iSession) = D4.cueT1(:,iC,iSession)*(-1); % *reflip data
%             D6.cueT2(:,iC,iSession) = D4.cueT2(:,iC,iSession)*(-1);
%             D6.cueN(:,iC,iSession) = D4.cueN(:,iC,iSession)*(-1);
%         else 
%             D6.cueT1(:,iC,iSession) = D4.cueT1(:,iC,iSession); 
%             D6.cueT2(:,iC,iSession) = D4.cueT2(:,iC,iSession);
%             D6.cueN(:,iC,iSession) = D4.cueN(:,iC,iSession);
%         end
%     end
% end

%% difference max amp 

% fields = {'cueT1','cueT2'}; 

figure
set(gcf,'Position',[100 100 1000 100])

for iT = 1:numel(targetNames)
    target = targetNames(iT);
    idxTargetName = find(strcmp(p.eventNames,target));
    idxT = find(t==p.eventTimes(idxTargetName));
    window = idxT:idxT+windowPostTarget+windowSize;
    for iC = 1:numel(channelGroups)
        channelGroup = sprintf('chTop%d',channelGroups(iC)); 
        
        [h, pval(iC), ci, stats] = ttest(groupA.(channelGroup).(targetNames{iT}).(fields{1}).max,groupA.(channelGroup).(targetNames{iT}).(fields{2}).max,'alpha',0.05);
        difference = groupA.(channelGroup).(targetNames{iT}).(fields{1}).max - groupA.(channelGroup).(targetNames{iT}).(fields{2}).max; 
        differenceMean = nanmean(difference);
        differenceSte = nanstd(difference)/sqrt(numel(subjectNames));
        
        subplot (1,numel(targetNames),iT)
        hold on
        % vertical error bar
        errV = errorbar(iC,...
            differenceMean,... % Y: max
            differenceMean-ci(1),ci(2)-differenceMean,... % vertical error bar, 95% confidence interval
            'Color',p.cueColors(4,:),'LineWidth',2,'CapSize',10);
        p1 = plot(iC,...
            differenceMean,...
            'o','MarkerSize',8);
        p1.Color = p.cueColors(4,:);
        p1.MarkerFaceColor = p.cueColors(4,:);
        xlim([0,numel(channelGroups)+1])
        xlabel('Number of channels')
        xticks(1:numel(channelGroups))
        xticklabels(string(num2cell(channelGroups)))
        ylim([-6.5e-14 6.5e-14])
        ylabel('Amp difference (T)')
        hline(0,':k')
        title(sprintf('%s max',targetNames{iT}))
    end
    disp(pval)
end
if saveFigs
    saveas(gcf,sprintf('%s/difference_maxAmp.svg',figDir)) 
end

%% difference max idx 

% fields = {'cueT1','cueT2'}; 

figure
set(gcf,'Position',[100 100 1000 100])

for iT = 1:numel(targetNames)
    target = targetNames(iT);
    idxTargetName = find(strcmp(p.eventNames,target));
    idxT = find(t==p.eventTimes(idxTargetName));
    window = idxT:idxT+windowPostTarget+windowSize;
    for iC = 1:numel(channelGroups)
        channelGroup = sprintf('chTop%d',channelGroups(iC)); 
        
        [h, pval(iC), ci, stats] = ttest(groupA.(channelGroup).(targetNames{iT}).(fields{1}).max,groupA.(channelGroup).(targetNames{iT}).(fields{2}).max,'alpha',0.05);
        difference = groupA.(channelGroup).(targetNames{iT}).(fields{1}).maxIdx - groupA.(channelGroup).(targetNames{iT}).(fields{2}).maxIdx; 
        differenceMean = nanmean(difference);
        differenceSte = nanstd(difference)/sqrt(numel(subjectNames));
        
        subplot (1,numel(targetNames),iT)
        hold on
        % vertical error bar
        errV = errorbar(iC,...
            differenceMean,... % Y: max
            differenceMean-ci(1),ci(2)-differenceMean,... % vertical error bar, 95% confidence interval
            'Color',p.cueColors(4,:),'LineWidth',2,'CapSize',10);
        p1 = plot(iC,...
            differenceMean,...
            'o','MarkerSize',8);
        p1.Color = p.cueColors(4,:);
        p1.MarkerFaceColor = p.cueColors(4,:);
        xlim([0,numel(channelGroups)+1])
        xlabel('Number of channels')
        xticks(1:numel(channelGroups))
        xticklabels(string(num2cell(channelGroups)))
        ylim([-20 20])
        ylabel('Latency difference (ms)')
        hline(0,':k')
        title(sprintf('%s',targetNames{iT}))
    end
    disp(pval)
end
if saveFigs
    saveas(gcf,sprintf('%s/difference_maxIdx.svg',figDir)) 
end

%% plot ERF 

figure
set(gcf,'Position',[100 100 1000 400])
hold on 
window = 1:4001; 
for iF = 1:numel(fields)
    plot(t(window),D.group.(fields{iF})(window),'LineWidth',2,'Color',p.cueColors(iF,:))
end
hline(0,':k')
vline(p.eventTimes,':k',p.eventNames) 
for iT = 1:numel(targetNames)
    target = targetNames(iT);
    idxTargetName = find(strcmp(p.eventNames,target));
    idxT = find(t==p.eventTimes(idxTargetName));
    window = idxT+windowPostTarget:idxT+windowSize; 
    l = line([t(window(1)),t(window(end))],[0,0]); 
    l.Color = [0 0 0 0.5]; 
    l.LineWidth = 5;
end
ylabel('Amplitude (T)')
xlabel('Time (ms)') 
legend(fields)
if saveFigs
    saveas(gcf,sprintf('%s/ERF.svg',figDir)) 
end


