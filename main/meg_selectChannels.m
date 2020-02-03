function C = meg_selectChannels(sessionDir,data)

% MEG_SELECTCHANNELS(sessionDir)
%
% INPUTS
% sessionDir 
% data 
%   matrix time x channel x trial 
% 
% OUPUT
%   C
%       structure of top channels by analysis type 
%
% Karen Tian
% January 2020

% snr, peak, topo, 

%% setup

% sessionDir = 'R0817_20190625';

exptShortName = 'TA2';
analStr = 'bietfp';
exptDir = '/Users/kantian/Dropbox/Data/TA2/MEG'; 

p = meg_params('TA2');

fileBase = sessionDirToFileBase(sessionDir, exptShortName);

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);

filename = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file* 
figDir = sprintf('%s/figures/%s', dataDir);

behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
behavFile = dir(sprintf('%s/*.mat', behavDir));
behav = load(sprintf('%s/%s', behavDir, behavFile.name));

eyesClosedBase = sessionDirToFileBase(sessionDir, 'EyesClosed');
eyesClosedFile = sprintf('%s/%s.sqd', dataDir, eyesClosedBase); 
eyesClosedFileBI = sprintf('%s/%s_bi.sqd', dataDir, eyesClosedBase); 

%% data info

sz = size(data); 
nT = sz(1); 
nCh = sz(2); 
nTrial = sz(3); 

t = p.tstart:p.tstop;
xlims = [min(t),max(t)]; % epoch time in ms

%% selection method 

selectPeak = 1; % sorts channels by ERF peak 
selectAlpha = 0; % sorts channels by alpha power 

%% select based off top peak 
meanVal = nanmean(data,3); % avg across trial
meanMeanVal = nanmean(meanVal); % avg across trial, then avg across time 

normVal = meanVal'./meanMeanVal'-1; 
colorbar
caxis([-100 100])

spacer = 5e-14; 

figure
for iC = peak.megChannels
    hold on
    plot(t, meanVal(:,iC) + spacer*(iC-1)) 
end
title('mean ch val')
xlim(xlims)
xlabel('time (ms)')
ylabel('amplitude')
vline(peak.eventTimes,'k',peak.eventNames)

for iC = peak.megChannels 
    val = meanVal(); 
end

 
figure
imagesc(meanVal')
colorbar

%     spacer = 40;
%     figure
%     for iF=1:nFields
%         vals = data.(fieldName{iF});
%         subplot (nFields,1,iF)
%         hold on
%         for iC = 1:nAllChannels
%             meanTrial = nanmean(vals(:,iC,:),3);
%             plot(t, abs(meanTrial) + spacer*iC)
%         end
%         title(sprintf('%s',fieldName{iF}))
%         xlim(xlims)
%         xlabel('time (ms)')
%         ylabel('amplitude')
%         vline(p.eventTimes,'k',p.eventNames)
%     end

%% peak selector

vals = nanmean(data,3);
nPk = 3;
twindow = p.eventTimes(2):p.eventTimes(3)+300; % window to look for peaks, from T1 to T2+300

for iC = 1:nCh
    if nanmean(vals(:,iC),1)==0 % why getting flat channels? check preproc bad channel output
        Pk.peaksPos(iC,:) = NaN;
        Pk.peaksPosVals(iC,:) = NaN;
        Pk.peaksNeg(iC,:) = NaN;
        Pk.peaksNegVals(iC,:) = NaN;
    else
        % positive peaks
        [pks,locs,~,prom] = findpeaks(vals(twindow,iC));
        [~, idx] = sort(prom,1,'descend');
        idx = sort(idx(1:nPk));
        peaksPos = twindow(locs(idx));
        peaksPosVals = vals(peaksPos,iC);
        
        % negative peaks
        [pks,locs,~,prom] = findpeaks(-vals(twindow,iC));
        [~, idx] = sort(prom,1,'descend');
        idx = sort(idx(1:nPk));
        peaksNeg = twindow(locs(idx));
        peaksNegVals = vals(peaksNeg,iC);
        
        Pk.peaksPos(iC,:) = peaksPos;
        Pk.peaksPosVals(iC,:) = peaksPosVals;
        Pk.peaksNeg(iC,:) = peaksNeg;
        Pk.peaksNegVals(iC,:) = peaksNegVals;  
    end
end

% subplot sizing 
nRow = 10; 
nBreak = 60; 
nCol = ceil(nBreak/nRow); 

figure
set(gcf,'Position',[100 100 2400 1200])
hold on
nFig = 1; 
for iC = 1:nCh
    if mod(iC,nBreak) == 0 % new figure at breakpoints 
        figure
        set(gcf,'Position',[100 100 2400 1200])
        hold on
        nFig = nFig + 1; 
    end
    subplot(nRow,nCol,iC-(nBreak*(nFig-1))+1)
    hold on
    plot(t,vals(:,iC))
    plot(Pk.peaksPos(iC,:), Pk.peaksPosVals(iC,:), '.g', 'MarkerSize', 30)
    plot(Pk.peaksNeg(iC,:), Pk.peaksNegVals(iC,:), '.b', 'MarkerSize', 30)
    vline(p.eventTimes,'k',p.eventNames)
    xlim(xlims)
    title(sprintf('%d',iC))
end


% figure
% plot(t(locs), p, '.') % visualize peak prominence



% store results of peaks analysis
peaks.measure = m;
peaks.t = t;
peaks.vals = vals;
peaks.nPk = nPk;
peaks.peaksPos = peaksPos;
peaks.peaksPosVals = peaksPosVals;
peaks.peaksNeg = peaksNeg;
peaks.peaksNegVals = peaksNegVals;


%% select based off top alpha channels 

if selectAlpha
    load(sprintf('%s/Alpha.mat',matDir));
    
    power = Alpha.TFRhann.powspctrm;
    powerAvg = nanmean(power,3); % average power across time
    alpha = nanmean(powerAvg(:,Alpha.tfAlpha-2:Alpha.tfAlpha+2),2); % individual alpha ± 2Hz, avg across time
    [C.powSortAlpha,C.chSortAlpha] = sort(alpha,'descend');
    
    figure
    set(gcf,'Position',[100 100 300 500])
    imagesc(alpha); % check
    set(gca,'xtick',[],'ytick',0:10:160)
    ylabel('channel')
    title(sprintf('%s eyes closed alpha %d ± 2Hz \n top 5 channels [%d %d %d %d %d]',und2space(sessionDir),Alpha.tfAlpha,C.chSortAlpha(1:5)))
end

%% save figs

fH = sort(findobj('Type','figure'));
figNames = {'TopChannels_Alpha'};
rd_saveAllFigs(fH, figNames, sessionDir, figDir, [])

%% save top channels 

save(sprintf('%s/C.mat',matDir),'C')

end



