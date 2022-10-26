function [A,fH,figNames] = meg_plotGroupITPC(data,sessionDir,p,selectedChannels,selectedFreq)

% MEG_PLOTGROUPITPC(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% Karen Tian
% July 2020

%% setup

condNames = {'all'}; % {'cueT1','cueT2'}; 
p = meg_params('TANoise_ITPC'); 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions('TANoise'); 

saveFigs = 1; 

figDir = '/Users/kantian/Dropbox/Data/TANoise/MEG/Group/figures/ITPCUpDown/phaseAngle';
if ~exist(figDir,'dir')
    mkdir(figDir)
end

%% params 

% timing 
t = p.tstart:p.tstop;
toi = t; 
xtick = 1:80:numel(toi);

% figures
foi = 20; 
ytick = 10:10:numel(foi);
xlims = [size(toi,1),size(toi,2)];

singleTrial = 0;
avgTrial = 1; 

%% make groupA by ITPC dir 

for i = 1:numel(sessionNames)
    for iC = 1:numel(condNames)
        if singleTrial
            groupITPC.(condNames{iC}).selectedChannels(i,:)= groupA{i}.(condNames{iC}).selectedChannels;
            groupITPC.(condNames{iC}).tfPows(:,i) = nanmean(groupA{i}.(condNames{iC}).tfPows,2);
            groupITPC.(condNames{iC}).tfPowsAvg(:,i) = nanmean(groupA{i}.(condNames{iC}).tfPowsAvg,2);
            groupITPC.(condNames{iC}).ITPC(:,i) = groupA{i}.(condNames{iC}).ITPCMean;
            groupITPC.(condNames{iC}).phaseAngle(i,:) = nanmean(groupA{i}.(condNames{iC}).phaseAngleCh,1);
            groupITPC.(condNames{iC}).phaseAngleCh(i,:,:) = groupA{i}.(condNames{iC}).phaseAngleCh;
            groupITPC.(condNames{iC}).phaseAngleTrial(i,:,:,:) = groupA{i}.(condNames{iC}).phaseAngle;
            
        elseif avgTrial
            groupITPC.(condNames{iC}).selectedChannels(i,:)= groupA{i}.(condNames{iC}).selectedChannels;
            groupITPC.(condNames{iC}).tfAmps(:,:,i) = groupA{i}.(condNames{iC}).tfAmps;
            groupITPC.(condNames{iC}).tfPows(:,:,i) = groupA{i}.(condNames{iC}).tfPows;
            groupITPC.(condNames{iC}).ITPC(:,:,i) = groupA{i}.(condNames{iC}).ITPC;
            groupITPC.(condNames{iC}).ITPCMean(:,i) = groupA{i}.(condNames{iC}).ITPCMean;
            groupITPC.(condNames{iC}).phaseAngle(:,:,i) = groupA{i}.(condNames{iC}).phaseAngle;
        end
    end
end

%% average to subjects 

nSubject = 1; 
for i = 1:2:numel(sessionNames)
    for iC = 1:numel(condNames)
        groupITPC.(condNames{iC}).tfPowsSubject(:,nSubject) = nanmean(groupITPC.(condNames{iC}).tfPows(:,i:i+1),2); 
        groupITPC.(condNames{iC}).tfPowsAvgSubject(:,nSubject) = nanmean(groupITPC.(condNames{iC}).tfPowsAvg(:,i:i+1),2); 
        groupITPC.(condNames{iC}).ITPCSubject(:,nSubject) = nanmean(groupITPC.(condNames{iC}).ITPC(:,i:i+1),2); 
        groupITPC.(condNames{iC}).phaseAngleSubject(nSubject,:) = nanmean(groupITPC.(condNames{iC}).phaseAngle(i:i+1,:),1); 
    end
    nSubject = nSubject + 1; 
end

%% plot avgtrial phase angle per session per 5 channels
% July 25, 2020

ITPCDirs = [1,-1];
colors = [97/255, 159/255, 145/255;...
    179/255, 117/255, 195/255];
timeOI = p.eventTimes(2)+abs(p.tstart); % T1 

figure
set(gcf, 'Position',  [100, 100, 1200, 500])
nPlot = 1; 
for iDir = 1:2
    ITPCoi = ITPCDirs(iDir);
    sessions = 1:numel(sessionNames);
    sessionIdx = sessions(ITPCsession==ITPCoi); % uppers 1, downers -1
    count = 1; 
    for iC = 1:numel(condNames)
        for iS = 1:numel(sessionIdx)
        for chOI = 1:5  
            theta(count) = groupITPC.(condNames{iC}).phaseAngle(timeOI,chOI,sessionIdx(iS));
            %             for i = 1:numel(sessionIdx)
            %                 thetaSession = [];
            %                 thetaSession = theta(i,:);
            %                 avgTheta(i) = circ_mean(thetaSession(~isnan(thetaSession)),[],2);
            %                 pp = polarplot([avgTheta(i),avgTheta(i)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',1);
            %                 pp.Color(4) = 0.1*i;
            %                 hold on
            %             end
            %     avgAngle = circ_mean(avgTheta,[],2);
            %     pp = polarplot([avgAngle,avgAngle],[0,1],'Color',p.cueColors(iC,:),'LineWidth',5);
            %     pp.Color(4) = 0.7;
            pp = polarplot([theta(count),theta(count)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',1);
            pp.Color = colors(iDir,:);
            pp.Color(4) = 0.5; 
            hold on 
            pflip = polarplot([theta(count)+pi,theta(count)+pi],[0,1],'Color',p.cueColors(iC,:),'LineWidth',1);
            pflip.Color = colors(iDir,:);
            pflip.Color(4) = 0.5; 
            % pause(0.5)
            count = count + 1; 
        end
        grid off
        ppAx = gca;
        ppAx.RColor = [0 0 0];
        ppAx.LineWidth = 1;
        % ppAx.ThetaTick = 0:30:330;
        % ppAx.ThetaTick = [];
        ppAx.RTick = [0,1];
        rlim([0,1])
        title(sprintf('Phase Angle at T1 \n per session per channel'))
        end
    end
end

% add average 
% avgAngle = circ_mean(testTheta(1,:),[],2);
% avgVec = circ_r(testTheta(1,:),[],[],2); 

if saveFigs 
    titleText = 'phaseAngle_T1'; 
    print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
end

%% average to sessions ITPC phase angle
% August 3, 2020

% subset and save T1 phase angle 
for i = 1:numel(sessions)
    phaseAngleT1 = squeeze(groupITPC.all.phaseAngle(timeOI,:,i));
    groupITPC.all.phaseAngleT1(:,i) = phaseAngleT1; 
end

%% test circ_mean

testThetaVec = 0:pi/4:2*pi-pi/2; 
testTheta = repmat(testThetaVec,2,1); 

testRho = [0,1];
testRho = repmat(testRho,size(testTheta,2),1); 
testRho = testRho'; 

avgAngle = circ_mean(testTheta(1,:),[],2);
avgVec = circ_r(testTheta(1,:),[],[],2); 

figure
pp = polarplot(testTheta,testRho,'Color',p.cueColors(1,:),'LineWidth',1);
hold on 
pp = polarplot([avgAngle,avgAngle],[0,avgVec],'Color',p.cueColors(2,:),'LineWidth',5);
titleText = 'Test circular mean (1)'; 
title(titleText)

if saveFigs 
    print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
end

%% save theta.Up/Down
% July 25, 2020 12pm 

ITPCDirs = [1,-1];
ITPCDirNames = {'Up','Down'}; 
% colors = [97/255, 159/255, 145/255;...
%     179/255, 117/255, 195/255];
timeOI = p.eventTimes(2)+abs(p.tstart); % T1 

figure
set(gcf, 'Position',  [100, 100, 1200, 500])
for iDir = 1:2
    ITPCoi = ITPCDirs(iDir);
    sessions = 1:numel(sessionNames);
    sessionIdx = sessions(ITPCsession==ITPCoi); % uppers 1, downers -1
    count = 1; 
    for iC = 1:numel(condNames)
        for iS = 1:numel(sessionIdx)
        for chOI = 1:5  
            theta.(ITPCDirNames{iDir})(count) = groupITPC.(condNames{iC}).phaseAngle(timeOI,chOI,sessionIdx(iS));
            pp = polarplot([theta.(ITPCDirNames{iDir})(count),theta.(ITPCDirNames{iDir})(count)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',1);
            pp.Color = colors(iDir,:);
            pp.Color(4) = 0.5; 
            hold on 
            count = count + 1; 
        end
        grid off
        ppAx = gca;
        ppAx.RColor = [0 0 0];
        ppAx.LineWidth = 1;
        % ppAx.ThetaTick = 0:30:330;
        % ppAx.ThetaTick = [];
        ppAx.RTick = [0,1];
        rlim([0,1])
        title(sprintf('Phase Angle at T1 \n per session per channel'))
        end
    end
end
if saveFigs 
    titleText = 'phaseAngle_T1_unflipped'; 
    print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
end

%% mirror phase angles, save to A 
% July 25, 2020 1pm 

colors = [0.4660, 0.6740, 0.1880, 0.5;...
     0.4940, 0.1840, 0.5560, 0.5];  

figure
for iDir = 1:numel(ITPCDirNames)
    thetaMirrorVec = theta.(ITPCDirNames{iDir}); 
    chs = size(thetaMirrorVec,2);
    for i = chs+1:chs*2
        thetaMirrorVec(:,i) = thetaMirrorVec(:,i-chs)+pi;
    end
    thetaMirror = repmat(thetaMirrorVec,2,1); 
    
    rho = [0,1]; 
    rho = repmat(rho,size(thetaMirrorVec,2),1); 
    rho = rho'; 
    
    
    pp = polarplot(thetaMirror,rho,'Color',colors(iDir,:),'LineWidth',1);
    hold on
    
    
    A.(ITPCDirNames{iDir}).thetaMirror = thetaMirror;
    A.(ITPCDirNames{iDir}).rho = rho; 

    title(sprintf('Phase Angle at T1 \n per session per channel'))
end
if saveFigs 
    titleText = 'phaseAngle_T1_mirror'; 
    print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
end

%% slice A back to 180 

for iDir = 1:numel(ITPCDirNames)
    vals = A.(ITPCDirNames{iDir}).thetaMirror;
    for i = 1:size(vals,2) % 80, or 100 
        for j = 1:size(vals,1) % 2
            if vals(j,i) < -pi/2 || vals(j,i) > pi/2
                vals(j,i) = NaN; 
            end
        end
    end
    newVals = []; 
    newVals = vals(~isnan(vals))';
    A.(ITPCDirNames{iDir}).thetaHalf = newVals; 
end

%% plot 180 
figure
for iDir = 1:numel(ITPCDirNames)
    angle = repmat(A.(ITPCDirNames{iDir}).thetaHalf,2,1); 
    
    rho = [0,1]; 
    rho = repmat(rho,size(angle,2),1); 
    rho = rho'; 
    
    % plot half 
    pp = polarplot(angle,rho,'Color',colors(iDir,:),'LineWidth',1);
    hold on
   
end

title(sprintf('Phase Angle at T1 \n per session per channel'))

if saveFigs 
    titleText = 'phaseAngle_T1_thetaHalf'; % 'phaseAngle_T1_doubleAvg'; 
    print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
end

%% double mirrored 180 to 360 

for iDir = 1:numel(ITPCDirNames)
    vals = []; 
    vals = A.(ITPCDirNames{iDir}).thetaHalf*2; 
    A.(ITPCDirNames{iDir}).thetaDouble = vals; 
end

%% doubled 
figure
for iDir = 1:numel(ITPCDirNames)
    angle = repmat(A.(ITPCDirNames{iDir}).thetaDouble,2,1); 
    
    rho = [0,1]; 
    rho = repmat(rho,size(angle,2),1); 
    rho = rho'; 
    
    % plot half 
    pp = polarplot(angle,rho,'Color',colors(iDir,:),'LineWidth',1);
    hold on
   
end

% calc mean 
avgUp= circ_mean(A.Up.thetaDouble,[],2);
vecUp = circ_r(A.Up.thetaDouble,[],[],2); 
confUp = circ_confmean(A.Up.thetaDouble');

avgDown= circ_mean(A.Down.thetaDouble,[],2);
vecDown = circ_r(A.Down.thetaDouble,[],[],2); 
confDown = circ_confmean(A.Down.thetaDouble');

% plot avg
% pp = polarplot([avgUp,avgUp],[0,vecUp],'Color',colors(1,:),'LineWidth',8);
% pp = polarplot([avgUp+confUp,avgUp+confUp],[0,1],'Color',colors(1,:),'LineWidth',1);
% pp = polarplot([avgUp-confUp,avgUp-confUp],[0,1],'Color',colors(1,:),'LineWidth',1);
% 
% pp = polarplot([avgDown,avgDown],[0,vecDown],'Color',colors(2,:),'LineWidth',8);
% pp = polarplot([avgDown+confDown,avgDown+confDown],[0,1],'Color',colors(2,:),'LineWidth',1);
% pp = polarplot([avgDown-confDown,avgDown-confDown],[0,1],'Color',colors(2,:),'LineWidth',1);

title(sprintf('Phase Angle at T1 \n per session per channel'))

if saveFigs 
    titleText = 'phaseAngle_T1_thetaDouble'; % 'phaseAngle_T1_doubleAvg'; 
    print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
end

%% circular stats 
% [corcc pval] = circ_corrcc(A.Up.thetaHalf, A.Down.thetaHalf(1:74)); % correlation? 

allAngles = [A.Up.thetaHalf A.Down.thetaHalf]; 
idxUp = 1:size(A.Up.thetaHalf,2); 
idxDown = size(A.Up.thetaHalf,2)+1:size(allAngles,2); 

% [pval table] = circ_hktest(allAngles, idxUp, idxDown, 0, ITPCDirNames); %
% needs equal lengths

% one way ANOVA for circular data 
[pval, table] = circ_wwtest(A.Up.thetaHalf, A.Down.thetaHalf);

%% avg plot + conf int

figure 
% up
pp = polarplot([avgUp,avgUp],[0,vecUp],'Color',colors(1,:),'LineWidth',8);
hold on 
pp = polarplot([avgUp+confUp,avgUp+confUp],[0,1],'Color',colors(1,:),'LineWidth',1);
pp = polarplot([avgUp-confUp,avgUp-confUp],[0,1],'Color',colors(1,:),'LineWidth',1);
% down 
pp = polarplot([avgDown,avgDown],[0,vecDown],'Color',colors(2,:),'LineWidth',8);
pp = polarplot([avgDown+confDown,avgDown+confDown],[0,1],'Color',colors(2,:),'LineWidth',1);
pp = polarplot([avgDown-confDown,avgDown-confDown],[0,1],'Color',colors(2,:),'LineWidth',1);

title(sprintf('Phase Angle at T1 \n per session per channel'))

if saveFigs 
    titleText = 'phaseAngle_T1_avgConfInt'; 
    print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
end

%% plot subjects single trial power 
figure
set(gcf, 'Position',  [100, 100, 800, 300])
hold on
for iC = 1:numel(condNames)
    errBar = shadedErrorBar(t,nanmean(groupITPC.(condNames{iC}).tfPowsSubject,2),nanstd(groupITPC.(condNames{iC}).tfPowsSubject,[],2)/sqrt(numel(subjectNames)));
    errBar.patch.FaceColor = p.cueColors(iC,:);
    errBar.edge(1).Color = p.cueColors(iC,:);
    errBar.edge(2).Color = p.cueColors(iC,:);
    plot(t,nanmean(groupITPC.(condNames{iC}).tfPowsSubject,2),'Color',p.cueColors(iC,:),'LineWidth',2)
end
vline(p.eventTimes,':k')
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
ylabel('Single Trial 20 Hz Power (fT^2)')
xlabel('Time (ms)')
title(sprintf('Group (subjects averaged)'))

%% plot subjects single trial power index upper downers 

ITPCoi = -1; 
subjects = 1:numel(subjectNames); 
whichSubjects = subjects(ITPCsubject==ITPCoi); % uppers 1, downers -1 

figure
set(gcf, 'Position',  [100, 100, 800, 300])
hold on
for iC = 1:numel(condNames)
    errBar = shadedErrorBar(t,nanmean(groupITPC.(condNames{iC}).tfPowsSubject(:,whichSubjects),2),nanstd(groupITPC.(condNames{iC}).tfPowsSubject(:,whichSubjects),[],2)/sqrt(numel(whichSubjects)));
    errBar.patch.FaceColor = p.cueColors(iC,:);
    errBar.edge(1).Color = p.cueColors(iC,:);
    errBar.edge(2).Color = p.cueColors(iC,:);
    plot(t,nanmean(groupITPC.(condNames{iC}).tfPowsSubject(:,whichSubjects),2),'Color',p.cueColors(iC,:),'LineWidth',2)
end
vline(p.eventTimes,':k')
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
ylabel('Single Trial 20 Hz Power (fT^2)')
xlabel('Time (ms)')
title(sprintf('Group (subjects ITPC Dir = %d averaged)',ITPCoi))

%% plot average trial power by upper downers 

ITPCoi = -1; 
subjects = 1:numel(subjectNames); 
whichSubjects = subjects(ITPCsubject==ITPCoi); % uppers 1, downers -1 

figure
set(gcf, 'Position',  [100, 100, 800, 300])
hold on
for iC = 1:numel(condNames)
    errBar = shadedErrorBar(t,nanmean(groupITPC.(condNames{iC}).tfPowsAvgSubject(:,whichSubjects),2),nanstd(groupITPC.(condNames{iC}).tfPowsAvgSubject(:,whichSubjects),[],2)/sqrt(numel(whichSubjects)));
    errBar.patch.FaceColor = p.cueColors(iC,:);
    errBar.edge(1).Color = p.cueColors(iC,:);
    errBar.edge(2).Color = p.cueColors(iC,:);
    plot(t,nanmean(groupITPC.(condNames{iC}).tfPowsAvgSubject(:,whichSubjects),2),'Color',p.cueColors(iC,:),'LineWidth',2)
end
vline(p.eventTimes,':k')
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
ylabel('Average Trial 20 Hz Power (fT^2)')
xlabel('Time (ms)')
title(sprintf('Group (subjects ITPC Dir = %d averaged)',ITPCoi))

%% plot ITPC 

ITPCoi = 1; 
subjects = 1:numel(subjectNames); 
whichSubjects = subjects(ITPCsubject==ITPCoi); % uppers 1, downers -1 

figure
set(gcf, 'Position',  [100, 100, 800, 300])
hold on
for iC = 1:numel(condNames)
    errBar = shadedErrorBar(t,nanmean(groupITPC.(condNames{iC}).ITPC(:,whichSubjects),2),nanstd(groupITPC.(condNames{iC}).ITPC(:,whichSubjects),[],2)/sqrt(numel(whichSubjects)));
    errBar.patch.FaceColor = p.cueColors(iC,:);
    errBar.edge(1).Color = p.cueColors(iC,:);
    errBar.edge(2).Color = p.cueColors(iC,:);
    plot(t,nanmean(groupITPC.(condNames{iC}).ITPC(:,whichSubjects),2),'Color',p.cueColors(iC,:),'LineWidth',2)
end
ylim([0 0.6])
vline(p.eventTimes,':k')
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
ylabel('ITPC')
xlabel('Time (ms)')
title(sprintf('Group (subjects ITPC Dir = %d averaged)',ITPCoi))

%% itpc polar plot angle 

ITPCoi = -1; 
subjects = 1:numel(subjectNames); 
whichSubjects = subjects(ITPCsubject==ITPCoi); % uppers 1, downers -1 

figure
set(gcf, 'Position',  [100, 100, 300, 400])
for iC = 1:numel(condNames)
    theta = groupITPC.(condNames{iC}).phaseAngleSubject(whichSubjects,p.eventTimes(2)+abs(p.tstart))';
    for i = 1:numel(whichSubjects)
        pp = polarplot([theta(i),theta(i)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',0.5);
        pp.Color(4) = 0.5; 
        hold on
    end
    pp = polarplot([circ_mean(theta,[],2),circ_mean(theta,[],2)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',5);
    pp.Color(4) = 0.7; 
end
grid off 
ppAx = gca;
ppAx.RColor = [0 0 0];
ppAx.LineWidth = 1; 
% ppAx.ThetaTick = 0:30:330; 
% ppAx.ThetaTick = []; 
ppAx.RTick = [0,1];
rlim([0,1])
title(sprintf('Phase Angle at T1 (ITPC Dir = %d)',ITPCoi))

%% itpc polar plot angle abs angle?? 

ITPCoi = 1; 
subjects = 1:numel(subjectNames); 
whichSubjects = subjects(ITPCsubject==ITPCoi); % uppers 1, downers -1 

figure
set(gcf, 'Position',  [100, 100, 300, 400])
for iC = 1:numel(condNames)
    theta = groupITPC.(condNames{iC}).phaseAngleSubject(whichSubjects,p.eventTimes(2))';
    for i = 1:numel(whichSubjects)
        pp = polarplot([theta(i),theta(i)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',0.5);
        pp.Color(4) = 0.5; 
        hold on
    end
    pp = polarplot([circ_mean(theta,[],2),nanmean(theta,2)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',5);
    pp.Color(4) = 0.7; 
end
grid off 
ppAx = gca;
ppAx.RColor = [0 0 0];
ppAx.LineWidth = 1; 
% ppAx.ThetaTick = 0:30:330; 
% ppAx.ThetaTick = []; 
ppAx.RTick = [0,1];
rlim([0,1])
title(sprintf('Phase Angle at T1 (ITPC Dir = %d)',ITPCoi))

%% plot subjects overlaid itpc phase angle

ITPCoi = -1; 
subjects = 1:numel(subjectNames); 
whichSubjects = subjects(ITPCsubject==ITPCoi); % uppers 1, downers -1 

thetaMirrorVec = groupITPC.(condNames{iC}).phaseAngleSubject(whichSubjects,:); 
toi = p.eventTimes(2)-100:p.eventTimes(3)+100; 
figure
set(gcf, 'Position',  [100, 100, 800, 300])
plot(toi,thetaMirrorVec(:,toi+abs(p.tstart)))
vline(p.eventTimes,':k')
xlim([toi(1),toi(end)])
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
ylabel('Phase angle (rad)')
xlabel('Time (ms)')
title(sprintf('subjects ITPC Dir = %d',ITPCoi))

%% plot subjects overlaid itpc

ITPCoi = 1; 
subjects = 1:numel(subjectNames); 
whichSubjects = subjects(ITPCsubject==ITPCoi); % uppers 1, downers -1 

thetaMirrorVec = groupITPC.(condNames{iC}).tfPowsAvgSubject(:,whichSubjects); 
toi = p.eventTimes(2)-100:p.eventTimes(3)+100; 
figure
set(gcf, 'Position',  [100, 100, 800, 300])
plot(toi,thetaMirrorVec(toi+abs(p.tstart),:)')
vline(p.eventTimes,':k')
xlim([toi(1),toi(end)])
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
% ylabel('ITPC')
ylabel('Average trial power (fT^2)')
xlabel('Time (ms)')
title(sprintf('subjects ITPC Dir = %d',ITPCoi))

%% plot try single session 5 channels itpc phase angle around targets 

figure
set(gcf, 'Position',  [100, 100, 1500, 400])
sessions = 1:20; 
uppers = sessions(ITPCsession==1);
downers = sessions(ITPCsession==-1); 
for whichSession = downers % 1:numel(sessionNames)
    for iCh = 1:5 % :5
        subplot (1,5,iCh)
        for iC = 1:numel(condNames)
            theta = groupA{whichSession}.(condNames{iC}).phaseAngle(:,iCh,p.eventTimes(2)+abs(p.tstart))';
%             for i = 1:192
%                 % pp = polarplot([theta(i),theta(i)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',0.5);
%                 hold on
%                 pp.Color(4) = 0.5;
%                 hold on
%             end
            pp = polarplot([circ_mean(theta(~isnan(theta)),[],2),circ_mean(theta(~isnan(theta)),[],2)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',5);
            hold on 
            pp.Color(4) = 0.7;
            % pause(1)
        end
        grid off
        ppAx = gca;
        ppAx.RColor = [0 0 0];
        ppAx.LineWidth = 1;
        % ppAx.ThetaTick = 0:30:330;
        % ppAx.ThetaTick = [];
        ppAx.RTick = [0,1];
        rlim([0,1])
        figName = sprintf('Phase Angle at T1 (channel %d)',iCh);
        title(figName)
        % saveas(gcf,sprintf('%s_Ch%d_PhaseAngle.png',sessionNames{whichSession},iCh))
        % close all
    end
end

%% plot sessions top 1 channel overlaid phase angle 

ITPCDirs = [1,-1];

figure
set(gcf, 'Position',  [100, 100, 1200, 500])
nPlot = 1; 
for iDir = 1:2
    % ITPCoi = 1; % up (1) or down (-1)
    ITPCoi = ITPCDirs(iDir);
    % chOI = 2; % channel to plot
    sessions = 1:numel(sessionNames);
    sessionIdx = sessions(ITPCsession==ITPCoi); % uppers 1, downers -1
    
    for chOI = 1:5
        subplot (2,5,nPlot)
        for iC = 1:numel(condNames)
            avgTheta = []; theta = [];
            theta = squeeze(groupITPC.(condNames{iC}).phaseAngleTrial(sessionIdx,:,chOI,p.eventTimes(2)+abs(p.tstart)));
            for i = 1:numel(sessionIdx)
                thetaSession = [];
                thetaSession = theta(i,:);
                avgTheta(i) = circ_mean(thetaSession(~isnan(thetaSession)),[],2);
                pp = polarplot([avgTheta(i),avgTheta(i)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',1);
                pp.Color(4) = 0.1*i;
                hold on
            end
            %     avgAngle = circ_mean(avgTheta,[],2);
            %     pp = polarplot([avgAngle,avgAngle],[0,1],'Color',p.cueColors(iC,:),'LineWidth',5);
            %     pp.Color(4) = 0.7;
        end
        grid off
        ppAx = gca;
        ppAx.RColor = [0 0 0];
        ppAx.LineWidth = 1;
        % ppAx.ThetaTick = 0:30:330;
        % ppAx.ThetaTick = [];
        ppAx.RTick = [0,1];
        rlim([0,1])
        title(sprintf('Phase Angle at T1 \n by session top channel %d \n (ITPC Dir = %d)',chOI,ITPCoi))
        nPlot = nPlot+1; 
    end
end

