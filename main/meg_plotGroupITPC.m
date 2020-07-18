function [A,fH,figNames] = meg_plotGroupITPC(data,sessionDir,p,selectedChannels,selectedFreq)

% MEG_PLOTGROUPITPC(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% Karen Tian
% July 2020

%% setup

condNames = {'cueT1','cueT2'}; 
p = meg_params('TANoise_ITPC'); 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions('TANoise'); 

%% params 

% timing 
t = p.tstart:p.tstop;
xtick = 1:80:numel(toi);

% figures
ytick = 10:10:numel(foi);
xlims = [size(toi,1),size(toi,2)]; 

%% make groupA by ITPC dir 

for i = 1:numel(sessionNames)
    for iC = 1:numel(condNames)
        groupITPC.(condNames{iC}).selectedChannels(i,:)= groupA{i}.(condNames{iC}).selectedChannels;
        groupITPC.(condNames{iC}).tfPows(:,i) = nanmean(groupA{i}.(condNames{iC}).tfPows,2);
        groupITPC.(condNames{iC}).tfPowsAvg(:,i) = nanmean(groupA{i}.(condNames{iC}).tfPowsAvg,2);
        groupITPC.(condNames{iC}).ITPC(:,i) = groupA{i}.(condNames{iC}).ITPCMean;
        groupITPC.(condNames{iC}).phaseAngle(i,:) = nanmean(groupA{i}.(condNames{iC}).phaseAngleCh,1);
        groupITPC.(condNames{iC}).phaseAngleCh(i,:,:) = groupA{i}.(condNames{iC}).phaseAngleCh;
        groupITPC.(condNames{iC}).phaseAngleTrial(i,:,:,:) = groupA{i}.(condNames{iC}).phaseAngle;
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

vals = groupITPC.(condNames{iC}).phaseAngleSubject(whichSubjects,:); 
toi = p.eventTimes(2)-100:p.eventTimes(3)+100; 
figure
set(gcf, 'Position',  [100, 100, 800, 300])
plot(toi,vals(:,toi+abs(p.tstart)))
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

vals = groupITPC.(condNames{iC}).tfPowsAvgSubject(:,whichSubjects); 
toi = p.eventTimes(2)-100:p.eventTimes(3)+100; 
figure
set(gcf, 'Position',  [100, 100, 800, 300])
plot(toi,vals(toi+abs(p.tstart),:)')
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

