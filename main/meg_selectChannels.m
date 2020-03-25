function [fH,figNames,Pk] = meg_selectChannels(sessionDir,data)

% MEG_SELECTCHANNELS(sessionDir)
%
% INPUTS
% sessionDir 
% data 
%   matrix time x channel x trial 
%
% Karen Tian
% January 2020


%% channel selection method 

selectPeak = 0; % sorts channels by ERF peak
selectThreshold = 1; % selects channels that cross a threshold max(abs(amp))
selectAlpha = 0; % sorts channels by alpha power 

%% setup

% sessionDir = 'R0817_20190625';
exptShortName = 'TA2';
analStr = 'bietfp';

exptDir = '/Users/kantian/Dropbox/Data/TA2/MEG'; 

fileBase = sessionDirToFileBase(sessionDir, exptShortName);

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);

filename = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file* 
figDir = sprintf('%s/figures/%s', dataDir);

% behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
% behavFile = dir(sprintf('%s/*.mat', behavDir));
% behav = load(sprintf('%s/%s', behavDir, behavFile.name));

% eyesClosedBase = sessionDirToFileBase(sessionDir, 'EyesClosed');
% eyesClosedFile = sprintf('%s/%s.sqd', dataDir, eyesClosedBase); 
% eyesClosedFileBI = sprintf('%s/%s_bi.sqd', dataDir, eyesClosedBase); 

%% data info

sz = size(data); 
nT = sz(1); 
nCh = sz(2); 
nTrial = sz(3); 

vals = nanmean(data,3); % avg across trial 

p = meg_params('TA2_Analysis'); 
t = p.tstart:p.tstop;
xlims = [min(t),max(t)]; % epoch time in ms

% threshold by max(abs(amp))
toi = p.eventTimes(2):p.eventTimes(3)+300; % times during which to check if max amp crossed with ref to epoch timing 
toiIdx = toi-p.tstart; % reference to 1 timing 
nChOI = 7; 
percentileCh = (nChOI/nCh)*100; 

%% peak selector

if selectPeak
    nPk = 7;
    % tWindow = p.eventTimes(2):p.eventTimes(3)+300; % window to look for peaks, from T1 to T2+300
    % valsWindow = vals(tWindow,:);
    for iC = 1:nCh
        if nanmean(vals(:,iC),1)==0 % why getting flat channels? check preproc bad channel output
            Pk.peaksPos(iC,:) = NaN;
            Pk.peaksPosVals(iC,:) = NaN;
            Pk.peaksNeg(iC,:) = NaN;
            Pk.peaksNegVals(iC,:) = NaN;
        else
            % positive peaks
            [pks,locs,~,prom] = findpeaks(vals(:,iC));
            [~, idx] = sort(prom,1,'descend'); % idx sorted by desc prominence
            % idx = find(prom>std(prom)*thresh); peak past threshold s
            idx = sort(idx(1:nPk)); % idx ascending of top nPk prom
            peaksPos = t(locs(idx));
            peaksPosVals = pks(idx);
            
            % negative peaks
            [pks,locs,~,prom] = findpeaks(-vals(:,iC));
            [~, idx] = sort(prom,1,'descend');
            idx = sort(idx(1:nPk));
            peaksNeg = t(locs(idx));
            peaksNegVals = -pks(idx);
            
            Pk.peaksPos(iC,:)       = peaksPos;
            Pk.peaksPosVals(iC,:)   = peaksPosVals;
            Pk.peaksNeg(iC,:)       = peaksNeg;
            Pk.peaksNegVals(iC,:)   = peaksNegVals;
        end
    end
    
    % subplot sizing
    nRow = 10;
    nBreak = 35;
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
        plot(t,vals(:,iC),'LineWidth',1)
        plot(Pk.peaksPos(iC,:), Pk.peaksPosVals(iC,:), '.g', 'MarkerSize', 10)
        plot(Pk.peaksNeg(iC,:), Pk.peaksNegVals(iC,:), '.b', 'MarkerSize', 10)
        xlim(xlims)
        vline(p.eventTimes,'k',p.eventNames)
        yline(0,'--');
        % vline(Pk.peaksPos(iC,:),'g')
        % vline(Pk.peaksNeg(iC,:),'b')
        title(sprintf('channel %d',iC))
        % ylim([min(min(vals)),max(max(vals))])
    end
    
    % store results of peaks analysis
    Pk.t         = t;
    Pk.vals      = vals;
    Pk.nPk       = nPk;
    
end

%% select based on max(abs(amp)) crossing threshhold

if selectThreshold 
    
    toiY = zeros(length(toi),1); 
    absVals = abs(vals); 
    toiVals = absVals(toiIdx,:); 
    [toiMax, Idx] = max(toiVals,[],1); 
    maxIdx = toi(Idx); 
    
    [toiMaxSort, toiMaxChSort] = sort(toiMax,'descend'); 

    thresholdPrctile = 100-percentileCh; % 96.8153% for 5 channels 
    thresholdVal = prctile(toiMax,thresholdPrctile); % cutoff val
    passCh = find(toiMax > thresholdVal); 
    
    
    % plot bar check max(abs(amp)) channels
    figure
    set(gcf,'Position',[100 100 800 300])
    hold on
    b = bar(toiMax,'FaceColor','flat'); 
    b.CData = repmat([0.8 0.8 0.8],[nCh,1]);
    b.CData(passCh,:) = repmat([0 1 0],[nChOI,1]); % highlight channels > threshold 
    hline(thresholdVal,'k',sprintf('threshold = %d',thresholdVal))
    title(sprintf('channels [%s] above threshold',num2str(passCh))) 
    ylabel('max(abs(amp))')
    xlabel('channel')
    
    % plot check group across time
    figure
    hold on 
    plot(t,absVals)
    plot(toi,toiY,'k','LineWidth',3) % highlight black time of interest for checking treshold
    xlim(xlims)
    ylim([0,1.1*max(toiMax)])
    vline(p.eventTimes,'k',p.eventNames)
    hline(thresholdVal,'k',sprintf('threshold = %d',thresholdVal))
    xlabel('time (ms)')
    ylabel('abs(amp)')
    
    % subplot sizing
    nRow = 10;
    nBreak = 35;
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
        subplot(nRow,nCol,iC-(nBreak*(nFig-1))+1) % one subplot per channel 
        hold on
        plot(t,absVals(:,iC)) 
        plot(toi,toiY,'k','LineWidth',3) 
        xlim(xlims)
        ylim([0,1.1*max(toiMax)])
        vline(p.eventTimes,'k',p.eventNames)
        yline(thresholdVal,'r');
        if ismember(iC,passCh)
            yline(thresholdVal,'g'); % if channel passes threshold, turn line green
        end
        title(sprintf('channel %d',iC))
        rd_supertitle(sprintf('channels [%s] above threshold = %d',num2str(passCh),thresholdVal)) 
    end
   
    % plot only channels that pass threshold, sorted by max 
    figure
    set(gcf,'Position',[100 100 2400 1200])
    hold on
    for iC = 1:nChOI
        subplot(nChOI,1,iC) % one subplot per channel
        hold on
        plot(t,absVals(:,toiMaxChSort(iC)))
        plot(toi,toiY,'k','LineWidth',3)
        xlim(xlims)
        ylim([0,1.1*max(toiMax)])
        vline(p.eventTimes,'k',p.eventNames)
        yline(thresholdVal,'g');
        title(sprintf('channel %s',num2str(toiMaxChSort(iC))))
        rd_supertitle2(sprintf('channels above threshold = %d',thresholdVal)) 
    end

    disp(passCh)
    
    % save threshold info 
    T.passCh = passCh; % channels that pass threshold 
    T.nChOI = nChOI; % number of channels that pass threshold 
    T.thresholdVal = thresholdVal; % threshold value 
    T.thresholdPrctile = thresholdPrctile;  % threshold percentile of channels
    T.toi = toi; % time of interest to look for max(abs(amp))
    T.toiMax = toiMax; % max(abs(amp)) values within toi 
    T.toiMaxSort = toiMaxSort; % sorted toiMax 
    T.toiMaxChSort = toiMaxChSort; % channel idx of sorted toiMax
    T.maxIdx = maxIdx; % idx of max w reference to epoch timing 
    
end
    
%% topo 

% figure
% set(gcf,'Position',[100 100 800 800])
% meg_multiplot(vals', [], [], passCh') % ERF all channels avg trial 

figure
set(gcf,'Position',[100 100 800 800])
meg_topoplot(toiMax', [], [], passCh') % avg trial, avg time, highlight channels 
title(sprintf('channels [%s] above threshold = %d',num2str(passCh),thresholdVal))

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

fH = sort(double(findobj(0,'Type','figure')));

if selectThreshold
    figDirThresh = sprintf('%sthreshold%d',figDir,nChOI); 
    figNames = {'histMaxAmp','ERF_ThreshAllChannels','ERF_Thresh1','ERF_Thresh2','ERF_Thresh3','ERF_Thresh4','ERF_Thresh5','ERF_passCh','topoPassCh'};
    if ~exist(figDirThresh,'dir')
        mkdir(figDirThresh)
    end
    rd_saveAllFigs(fH, figNames, sessionDir, figDirThresh, [])
end

% figNames = {'TopChannels_Alpha'};
% figNames = {'ERFPeaks1','ERFPeaks2','ERFPeaks3','ERFPeaks4','ERFPeaks5','TopoERF'}
% rd_saveAllFigs(fH, figNames, sessionDir, figDir, [])

%% save info 

% save(sprintf('%s/Pk.mat',matDir),'Pk')
% save(sprintf('%s/C.mat',matDir),'C')

if selectThreshold
    save(sprintf('%s/T.mat',matDir),'T') % top channels by thresholding
end

%% close 
close all

end



