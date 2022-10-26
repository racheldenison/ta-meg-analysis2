% function [fH,figNames] = meg_plotITPCTimeSeries(sessionDir,i)

%% input 
saveFigs = 1;
user = 'karen'; 

expt = 'TANoise'; 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 

%% loop sessions 
for iS = 1 %:20
    sessionDir = sessionNames{iS};
    
    %% setup
    ch = groupC(iS).channelsRanked(1:5);
    p = meg_params('TANoise_ITPC');
    
    exptDir = meg_pathToTAMEG(expt, user);
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    
    ITPCmat = sprintf('%s/ITPC.mat',matDir);
    load(ITPCmat)
    t = 1:size(A.all.ITPC,1);
    
    analStr = 'ebi';
    figDir = sprintf('%s/figures/%s/ITPC5', dataDir, analStr);
    if ~exist(figDir,'dir')
        mkdir(figDir)
    end
    
    %% calculate ITPC by different windows
    
    vals = nanmean(A.all.ITPC(:,ch),2); % top 5 channels
    meanITPC = nanmean(A.all.ITPC,1); % avg ITPC across time per channel
    normITPC = A.all.ITPC - meanITPC; % norm ITPC by channel
    
    % pre precue, pre T1, target
    toi1 = -250+abs(p.tstart):0+abs(p.tstart); % preCue times
    toi2 = 750+abs(p.tstart):1000+abs(p.tstart); % preTarget
    toi3 = p.eventTimes(2)+abs(p.tstart):p.eventTimes(2)+abs(p.tstart)+250;
    
    vals1 = A.all.ITPC(toi1,:); % pre precue
    vals2 = A.all.ITPC(toi2,:); % pre T1
    vals3 = A.all.ITPC(toi3,:); % target
    
    normVals1 = normITPC(toi1,:);
    normVals2 = normITPC(toi2,:);
    normVals3 = normITPC(toi3,:);
    
    %% ITPC time series
    figure
    set(gcf,'Position',[100 100 800 300])
    vals_downsample = movmean(vals,10);
    % plot(t(1:10:end),vals(1:10:end),'k','LineWidth',2)
    hold on
    for i = 1:numel(p.eventTimes)
        xline(p.eventTimes(i)+abs(p.tstart),'Color',[0.5 0.5 0.5],'LineWidth',1)
    end
    
    % plot([toi1(1),toi1(end)],[0,0],'b','LineWidth',3)
    % plot([toi2(1),toi2(end)],[0,0],'g','LineWidth',3)
    % plot([toi3(1),toi3(end)],[0,0],'r','LineWidth',3)
    
    set(gca,'XTick',t(1):1000:t(end));
    set(gca,'XTickLabel',t(1)+p.tstart-1:1000:t(end)+p.tstart-1);
    
    plot(t,A.all.ITPC,'Color',[0.5 0.5 0.5 0.5],'LineWidth',0.2);
    p1 = plot(t(toi1),vals1,'Color',[0 0 1 0.3],'LineWidth',0.2);
    p2 = plot(t(toi2),vals2,'Color',[0 1 0 0.3],'LineWidth',0.2);
    p3 = plot(t(toi3),vals3,'Color',[1 0 0 0.3],'LineWidth',0.2);
    
    plot(t(1:10:end),vals_downsample(1:10:end),'k','LineWidth',3) % plot top 5 20 Hz channels
    ylim([0 0.8])
    xlabel('Time (ms)')
    ylabel('ITPC')
    title(sprintf('%s',und2space(sessionDir)))
    set(gca,'TickDir','out');
    ax = gca;
    ax.LineWidth = 1.5;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.FontSize = 14;
    box off
    set(gca,'TitleFontSizeMultiplier',1.3)
    
    %% topo fixed colorbar
    % topo pre preCue
    figure
    meg_topoplot(vals1',[],[],ch)
    title(sprintf('ITPC %s %d to %d ms',und2space(sessionDir),toi1(1)-abs(p.tstart),toi1(end)-abs(p.tstart)))
    caxis([0 0.5])
    colorbar
    
    % topo pre T1
    figure
    meg_topoplot(vals2',[],[],ch)
    title(sprintf('ITPC %s %d to %d ms',und2space(sessionDir),toi2(1)-abs(p.tstart),toi2(end)-abs(p.tstart)))
    caxis([0 0.5])
    colorbar
    
    % topo target times
    figure
    meg_topoplot(vals3',[],[],ch)
    title(sprintf('ITPC %s %d to %d ms',und2space(sessionDir),toi3(1)-abs(p.tstart),toi3(end)-abs(p.tstart)))
    caxis([0 0.5])
    colorbar
    
    %% variable colorbar
    % topo pre preCue
    figure
    meg_topoplot(vals1',[],[],ch)
    title(sprintf('ITPC %s %d to %d ms',und2space(sessionDir),toi1(1)-abs(p.tstart),toi1(end)-abs(p.tstart)))
    colorbar
    
    % topo pre T1
    figure
    meg_topoplot(vals2',[],[],ch)
    title(sprintf('ITPC %s %d to %d ms',und2space(sessionDir),toi2(1)-abs(p.tstart),toi2(end)-abs(p.tstart)))
    colorbar
    
    % topo target times
    figure
    meg_topoplot(vals3',[],[],ch)
    title(sprintf('ITPC %s %d to %d ms',und2space(sessionDir),toi3(1)-abs(p.tstart),toi3(end)-abs(p.tstart)))
    colorbar
    
    %% topo normalized by channel
    % topo pre preCue
    figure
    meg_topoplot(normVals1',[],[],ch)
    title(sprintf('normalized ITPC %s %d to %d ms',und2space(sessionDir),toi1(1)-abs(p.tstart),toi1(end)-abs(p.tstart)))
    caxis([-0.2 0.2])
    colorbar
    
    % norm topo pre T1
    figure
    meg_topoplot(normVals2',[],[],ch)
    title(sprintf('ITPC %s %d to %d ms',und2space(sessionDir),toi2(1)-abs(p.tstart),toi2(end)-abs(p.tstart)))
    caxis([-0.2 0.2])
    colorbar
    
    % norm topo target
    figure
    meg_topoplot(normVals3',[],[],ch)
    title(sprintf('ITPC %s %d to %d ms',und2space(sessionDir),toi3(1)-abs(p.tstart),toi3(end)-abs(p.tstart)))
    caxis([-0.2 0.2])
    colorbar
    
    %% return figure handle
    
    fH = sort(double(findobj(0,'Type','figure')));
    figNames = {'1_ITPCTimeSeries','2_ITPC_topo_preCue_fixedCAxis','3_ITPC_topo_preT1_fixedCAxis','4_ITPC_topo_T1_fixedCAxis',...
        '5_ITPC_topo_preCue_variableCAxis','6_ITPC_topo_preT1_variableCAxis','7_ITPC_topo_preCue_variableCAxis',...
        '8_ITPC_topo_preCue_norm','9_ITPC_topo_preT1_norm','10_ITPC_topo_preCue_norm'};
    if saveFigs
        rd_saveAllFigs(fH, figNames, [], figDir)
    end
    close all
end



