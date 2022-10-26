% meg_plotTopo

% MEG_PLOTTOPO
%
% Karen Tian
% August 2020

%% input 
expt = 'TANoise';
saveFigs = 1; 
figDir = '/Users/kantian/Dropbox/Data/TANoise/MEG/Group/figures/ITPCUpDown'; 
user = 'karen'; 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 

%% load data
% average data 
avgData = nanmean(data,3); 
D.all = avgData; 

%% setup 
p = meg_params(sprintf('%s_Analysis',expt)); 

%% plot 20Hz pows 
p = meg_params('TANoise_ITPC'); 
parentFigDir = '/Users/kantian/Dropbox/Data/TANoise/MEG/Group/figures/topo20HzPowsEps';
for iS = 1:numel(sessionNames)
    chs = groupC(iS).channelsRanked;
    chsHighlight = chs(1:5); 
    figDir = sprintf('%s/%s',parentFigDir,sessionNames{iS});
    if ~exist(figDir,'dir')
        mkdir(figDir)
    end
    figure
    count = 0;
    for iT = -250+abs(p.tstart):250:abs(p.tstart)+2000
        allVals = groupA{iS}.tfPows;
        maxVal = max(max(allVals));
        vals = groupA{iS}.tfPows(iT,:);
        meg_topoplot(vals',[],[],chsHighlight)
        caxis([0 maxVal*0.6])
        colorbar
        colormap('parula')
        title(sprintf('%s, %d ms',und2space(sessionNames{iS}),iT-abs(p.tstart)))
        titleText = sprintf('topo_%d_%s_%d',count,sessionNames{iS},iT-abs(p.tstart));
        print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
        count = count+1; 
    end
    close all
%     figure
%     cb = colorbar;
%     caxis([0 maxVal*0.6])
%     set(get(cb,'label'),'string','20 Hz Power');
%     titleText = sprintf('topo_%s_colorbar',sessionNames{iS});
%     print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
end

%% plot amp 
parentFigDir = '/Users/kantian/Dropbox/Data/TANoise/MEG/Group/figures/topoRawAmp';
for iS = 1:numel(sessionNames)
    chs = groupC(iS).channelsRanked;
    chsHighlight = chs(1:5); 
    figDir = sprintf('%s/%s',parentFigDir,sessionNames{iS});
    if ~exist(figDir,'dir')
        mkdir(figDir)
    end
    for iT = -250+500:250:500+2000
        figure
        vals = groupD{iS}.all(iT,:);
        meg_topoplot(vals',[],[],chsHighlight)
        colorbar
        title(sprintf('%s, %d ms',und2space(sessionNames{iS}),iT))
        titleText = sprintf('topo_%s_%d',sessionNames{iS},iT);
        print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
        close all
    end
end

%%
sessionIdx = 1; 
sessionDir = sessionNames{sessionIdx};

% time freq 
selectedChannels = p.megChannels; 
p = meg_params('TANoise_ITPC');

[A,fH,figNames] = meg_plotITPCavg(D,sessionDir,p,selectedChannels,20);
% [A,fH,figNames] = meg_plotTF(D,p,selectedChannels);

%% August 13, 2020
% plot average from -250 to 0ms of 20Hz power, ITPC, alpha  

variable = 'ITPC'; % 'tfPpws'; 'ITPC'
variableText = 'ITPC'; % '20 Hz Power'; 'ITPC'
parentFigDir = '/Users/kantian/Dropbox/Data/TANoise/MEG/Group/figures/ITPCUpDown/topo_20Hz_157Ch';
if ~exist(parentFigDir,'dir')
    mkdir(parentFigDir)
end
figDir = parentFigDir; 

% toi = -250+abs(p.tstart):0+abs(p.tstart); % preCue times 
toi = 750+abs(p.tstart):1000+abs(p.tstart); % preTarget 
% toi = p.eventTimes(2)+abs(p.tstart):p.eventTimes(2)+abs(p.tstart)+250; %
% target times 

%% loop subjects and save figs 
for iS = 1:2:numel(sessionNames)
    chs = groupC(iS).channelsRanked;
    chsHighlight = chs(1:5); % top 5 20Hz channels
    %     figDir = sprintf('%s/%s',parentFigDir,sessionNames{iS});
    %     if ~exist(figDir,'dir')
    %         mkdir(figDir)
    %     end
    vals = []; 
    vals = nanmean(groupA(iS).(variable)(toi,:),1); % average across time x ch 
    
    vals2 = nanmean(groupA(iS+1).(variable)(toi,:),1); % session 2
    valsBothSessions = nanmean(cat(1,vals,vals2),1); 
    vals = valsBothSessions; 
    
    figure
    meg_topoplot(vals',[],[],chsHighlight)
    colorbar
    % caxis([0,9e-13]) % amp 
    % caxis([0 15e-25]) % power 
    % caxis([0 0.5]) % ITPC
    
    title(sprintf('%s, %s, %d:%d ms',und2space(sessionNames{iS}),variableText,toi(1)-abs(p.tstart),toi(end)-abs(p.tstart)))
    titleText = sprintf('topo_%s_%d-%dms_%s',sessionNames{iS},toi(1)-abs(p.tstart),toi(end)-abs(p.tstart),variableText);
    set(gcf, 'Color', 'None')
    set(gca,'Color','none')
    % print(sprintf('%s/%s.eps',figDir,titleText),'-depsc') % eps format
    if saveFigs 
        print(sprintf('%s/%s.png',figDir,titleText),'-dpng')
    end
    close all
end









